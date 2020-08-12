import os
import numpy as np
import healpy as hp
import healsparse
from multiprocessing import Pool

import lsst.log

from .patch_mapper import PatchMapper, pool_initializer
from .utils import approx_patch_polygon_area, op_str_to_code


class MultiMapper(object):
    """
    Map a combination of tracts/filters/patches.
    """
    def __init__(self, butler, config, outputpath, ncores=1):
        """
        Instantiate a MultiMapper

        Parameters
        ----------
        butler : `lsst.daf.persistence.Butler`
           Data butler to retrieve exposures
        config : `Configuration`
           supreme configuration object
        outputpath : `str`
           Base path for output files
        ncores : `int`
           Number of cores to run with
        """
        self.butler = butler
        self.config = config
        self.outputpath = outputpath
        self.ncores = ncores

        lsst.log.setLevel("", lsst.log.ERROR)

        if not os.path.isdir(outputpath):
            raise RuntimeError("Outputpath %s does not exist." % (outputpath))

    def __call__(self, tracts, filters, patches=None, find_only=False, consolidate=True,
                 clobber=False, do_raise=False):
        """
        Compute the maps for a combination of tracts, filters, and patches.

        Will look for all combinations of tracts/filters/patches to generate
        maps.  If patches is None, will look for all possible patches for
        each tract/filter combination.

        Parameters
        ----------
        tracts : `list` of `int`
           List of tracts to look for coadds to generate maps
        filters : `list` of `str`
           List of filters to look for coadds to generate maps
        patches : `list` of `str`, optional
           List of patches to look for coadds to generate maps
        find_only : `bool`, optional
           Only find the data to run on.  Used for testing.
        consolidate : `bool`, optional
           Consolidate all patches for a given tract/filter into
           patch/filter maps for saving
        clobber : `bool`, optional
           Clobber any existing files
        do_raise : `bool`, optional
           Raise if there are any failures.  Otherwise continue to
           next tract / filter
        """
        skymap = self.butler.get('deepCoadd_skyMap')

        # We need to know the patch names
        if patches is None:
            patches = self._get_all_patch_names(skymap[tracts[0]])

        # Figure out the tract/filter/patch combinations
        multi_dict = {}
        for tract in tracts:
            for f in filters:
                for patch in patches:
                    if self.butler.datasetExists('deepCoadd', tract=tract,
                                                 filter=f, patch=patch):
                        if tract not in multi_dict:
                            multi_dict[tract] = {}
                        if f not in multi_dict[tract]:
                            multi_dict[tract][f] = []
                        multi_dict[tract][f].append(patch)

        for tract in multi_dict:
            for f in multi_dict[tract]:
                print('Found %d patches for tract %d / filter %s' %
                      (len(multi_dict[tract][f]), tract, f))

        if find_only:
            return multi_dict

        if len(multi_dict) == 0:
            print('No patch coadds found for any tracts / filters')
            return

        # Figure out the optimal nside_coverage per tract
        nside_coverage_tract = self._compute_nside_coverage_tract(skymap[tracts[0]])

        # Get the patch_mapper ready
        patch_mapper = PatchMapper(self.outputpath, self.config)

        # For each tract/filter we map to run all the patches
        for tract in multi_dict:
            for f in multi_dict[tract]:

                # Check for existing files, if clobber is False
                map_run = None
                ntot = 0
                nskip = 0
                if not clobber:
                    # This logic here only works for tract maps, and not for patch
                    # maps.  To support what I think is an edge case (where it
                    # would speed up skipping in some instances) is not worth it.
                    map_run = {}
                    for i, map_type in enumerate(self.config.map_types):
                        map_run[map_type] = []
                        for j, op_str in enumerate(self.config.map_types[map_type]):
                            ntot += 1
                            op_code = op_str_to_code(op_str)
                            fname = os.path.join(self.outputpath,
                                                 self.config.tract_relpath(tract),
                                                 self.config.tract_map_filename(f,
                                                                                tract,
                                                                                map_type,
                                                                                op_code))
                            map_run[map_type].append(not os.path.isfile(fname))
                            if not map_run[map_type][-1]:
                                nskip += 1
                                print('Found %s; skipping render as clobber is False' %
                                      (fname))

                    if nskip == ntot:
                        print('All maps have already been rendered and clobber is False.')
                        continue

                print('Running on tract %d / filter %s with %d cores.' %
                      (tract, f, self.ncores))

                values = zip([tract]*len(multi_dict[tract][f]),
                             [f]*len(multi_dict[tract][f]),
                             multi_dict[tract][f],
                             [consolidate]*len(multi_dict[tract][f]),
                             [map_run]*len(multi_dict[tract][f]),
                             [clobber]*len(multi_dict[tract][f]))

                pool = Pool(processes=self.ncores,
                            initializer=pool_initializer,
                            initargs=(self.butler, self.ncores == 1))
                results = pool.starmap(patch_mapper, values, chunksize=1)
                pool.close()
                pool.join()

                if not consolidate:
                    continue

                # Check for failures, raise if requested
                if np.any(np.array([r[2] for r in results])):
                    print("FAIL: Map failure on %d / %s" % (tract, f))
                    if do_raise:
                        raise RuntimeError("Exiting after failure on %d / %s" %
                                           (tract, f))
                    else:
                        continue

                for i, map_type in enumerate(self.config.map_types):
                    for j, op_str in enumerate(self.config.map_types[map_type]):
                        if map_run is not None:
                            if not map_run[map_type][j]:
                                continue

                        print('Consolidating tract %s / filter %s: %s, %s' %
                              (tract, f, map_type, op_str))

                        self._consolidate_tract_from_results(tract,
                                                             f,
                                                             nside_coverage_tract,
                                                             map_type,
                                                             op_str_to_code(op_str),
                                                             results, i, j,
                                                             clobber=clobber)

    def _consolidate_tract_from_results(self, tract, filtername, nside_coverage_tract,
                                        map_type, op_code,
                                        results, map_index, op_index,
                                        clobber=False):
        """
        Consolidate tract from the results.  Suitable to be run asynchronously.

        Parameters
        ----------
        tract : `int`
           Tract number
        filtername : `str`
           Name of the filter
        nside_coverage_tract : `int`
           nside_coverage for the tract map
        map_type : `str`
           Name of map type
        op_code : `int`
           Code for the operation
        results : `list` of `tuple`
           Parallelized results list
        map_index : `int`
           Index of map
        op_index : `int`
           Index of operation
        """
        dt = results[0][1][map_index][:, op_index].dtype

        tract_map = healsparse.HealSparseMap.make_empty(nside_coverage=nside_coverage_tract,
                                                        nside_sparse=self.config.nside,
                                                        dtype=dt)

        for patch_input_map, map_values_list, fails in results:
            tract_map.update_values_pix(patch_input_map.valid_pixels,
                                        map_values_list[map_index][:, op_index])

        fname = os.path.join(self.outputpath,
                             self.config.tract_relpath(tract),
                             self.config.tract_map_filename(filtername,
                                                            tract,
                                                            map_type,
                                                            op_code))
        tract_map.write(fname, clobber=clobber)

    def _compute_nside_coverage_tract(self, tract_info):
        """
        Compute the optimal coverage nside for a tract.

        Parameters
        ----------
        tract_info : `lsst.skymap.ExplicitTractInfo`
           Information on the tract

        Returns
        -------
        nside_coverage_patch : `int`
           Optimal coverage nside
        """
        num_patches = tract_info.getNumPatches()

        patch_area = approx_patch_polygon_area(tract_info.getPatchInfo((0, 0)),
                                               tract_info.getWcs())
        tract_area = patch_area * num_patches[0] * num_patches[1]
        nside_coverage_tract = 32
        while hp.nside2pixarea(nside_coverage_tract, degrees=True) > tract_area:
            nside_coverage_tract = int(2*nside_coverage_tract)
        nside_coverage_tract = int(np.clip(nside_coverage_tract / 2, 32, None))

        return nside_coverage_tract

    def _get_all_patch_names(self, tract_info):
        """
        Get all possible patch names for a tract.

        Parameters
        ----------
        tract_info : `lsst.skymap.ExplicitTractInfo`
           Information on the tract

        Returns
        -------
        patch_names : `list` of `str`
        """
        num_patches = tract_info.getNumPatches()

        patch_names = []
        for ii in range(num_patches[0]):
            for jj in range(num_patches[1]):
                patch_names.append('%d,%d' % (ii, jj))

        return patch_names
