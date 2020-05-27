import os
import numpy as np
import healpy as hp
import healsparse
# import multiprocessing
# from multiprocessing import Pool

from .patch_mapper import PatchMapper
from .utils import approx_patch_polygon_area, op_str_to_code


class MultiMapper(object):
    """
    Map a combination of tracts/filters/patches.
    """
    def __init__(self, butler, config, outputpath, ncores=1):
        """
        """
        self.butler = butler
        self.config = config
        self.outputpath = outputpath
        self.ncores = ncores

    def run(self, tracts, filters, patches=None, find_only=False):
        """
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
        patch_mapper = PatchMapper(self.butler, self.config, self.outputpath)

        # For each tract/filter we map to run all the patches
        for tract in multi_dict:
            for f in multi_dict[tract]:
                print('Running on tract %d / filter %s with %d cores.' %
                      (tract, f, self.ncores))

                # values = zip([f]*len(multi_dict[tract][f]),
                #              [tract]*len(multi_dict[tract][f]),
                #              multi_dict[tract][f])

                # pool = Pool(processes=self.ncores)
                # results = pool.starmap(patch_mapper.run, values, chunksize=1)
                # pool.close()
                # pool.join()
                results = [patch_mapper.run(f, tract, p) for p in multi_dict[tract][f]]

                for i, map_type in enumerate(self.config.map_types):
                    for j, op_str in enumerate(self.config.map_types[map_type]):
                        op_code = op_str_to_code(op_str)
                        patch_input_map = results[0][0]
                        map_values_list = results[0][1]

                        nct = nside_coverage_tract
                        ns = self.config.nside
                        dt = map_values_list[i][:, j].dtype
                        tract_map = healsparse.HealSparseMap.make_empty(nside_coverage=nct,
                                                                        nside_sparse=ns,
                                                                        dtype=dt)

                        # Put together all the patch maps
                        for patch_input_map, map_values_list in results:
                            tract_map.update_values_pix(patch_input_map.valid_pixels,
                                                        map_values_list[i][:, j])

                        fname = os.path.join(self.outputpath,
                                             self.config.tract_map_filename(f,
                                                                            tract,
                                                                            map_type,
                                                                            op_code))
                        tract_map.write(fname)

    def _compute_nside_coverage_tract(self, tract_info):
        """
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
        num_patches = tract_info.getNumPatches()

        patch_names = []
        for ii in range(num_patches[0]):
            for jj in range(num_patches[1]):
                patch_names.append('%d,%d' % (ii, jj))

        return patch_names
