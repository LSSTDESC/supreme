import os
import healsparse
import re

from .utils import op_str_to_code


class TractConsolidator(object):
    """
    Consolidate tract maps for a list of filters/maps.
    """
    def __init__(self, config, outputpath):
        """
        Instantiate a TractConsolidator.

        Parameters
        ----------
        config : `Configuration`
           supreme configuration object
        outputpath : `str`
           Base path for output files
        """
        self.config = config
        self.outputpath = outputpath

        if not os.path.isdir(outputpath):
            raise RuntimeError("Outputpath %s does not exist." % (outputpath))

    def __call__(self, filters, tracts=None, clobber=False, nside_coverage=None,
                 outputbase=None):
        """
        Consolidate maps for a combination of tracts / filters.

        Will use the configuration object to determine which maps need
        consolidation.

        Parameters
        ----------
        filters : `list` of `str`
           List of filters to look for maps to consolidate
        tracts : `list` of `int`, optional
           List of tracts to look for maps to consolidate.
           If None, will consolidate all tracts in the output directory.
        clobber : `bool`, optional
           Clobber any existing files
        nside_coverage : `int`, optional
           nside for the healsparse coverage map.  If None, use the same
           as the tracts.
        outputbase : `str`, optional
           Filename base for output consolidate files.

        Returns
        -------
        tracts : `list` of `int`
           List of tracts that were run.
        map_files : `list` of `str`
           List of map files that were created.
        map_inputs : `list` of `list`
           List of list of input files to each map.
        """
        if tracts is None:
            # Must find all the tracts available
            # These are the subdirectories that are numeric
            subdirs = [os.path.basename(d) for d in os.listdir(self.outputpath)
                       if os.path.isdir(os.path.join(self.outputpath, d)) and
                       re.search(r'^\d+$', d) is not None]
            if len(subdirs) == 0:
                print('No tract directories found in %s.' % (self.outputpath))
                return [], [], []
            tracts = [int(d) for d in subdirs]

        if outputbase is None:
            outputbase = self.config.outbase

        map_files = []
        map_inputs = []

        print('Consolidating maps from %d tracts' % (len(tracts)))
        for f in filters:
            for i, map_type in enumerate(self.config.map_types):
                for j, op_str in enumerate(self.config.map_types[map_type]):
                    op_code = op_str_to_code(op_str)

                    # Look for appropriate map files
                    consolidate_files = []
                    for tract in tracts:
                        fname = os.path.join(self.outputpath,
                                             self.config.tract_relpath(tract),
                                             self.config.tract_map_filename(f,
                                                                            tract,
                                                                            map_type,
                                                                            op_code))
                        if os.path.isfile(fname):
                            consolidate_files.append(fname)
                        else:
                            print("Warning: Tract %d does not have map %s" %
                                  (tract, os.path.basename(fname)))

                    if len(consolidate_files) == 0:
                        print("Could not find any tracts with map %s_%s." %
                              (map_type, op_str))
                        continue

                    print("Found %d tracts with map %s_%s to consolidate." %
                          (len(consolidate_files), map_type, op_str))

                    outfile = os.path.join(self.outputpath,
                                           self.config.consolidated_map_filename(outputbase,
                                                                                 f,
                                                                                 map_type,
                                                                                 op_code))
                    if os.path.isfile(outfile) and not clobber:
                        print("%s already exists and clobber is False.  Skipping..." %
                              (outfile))
                        continue

                    print("Consolidating to %s" % (outfile))

                    # Use in_memory mode because of compression/appending...
                    # Because we need to load the whole thing in memory anyway.
                    healsparse.cat_healsparse_files(consolidate_files, outfile,
                                                    check_overlap=False, clobber=clobber,
                                                    in_memory=True,
                                                    nside_coverage_out=nside_coverage)
                    map_files.append(outfile)
                    map_inputs.append(consolidate_files)

        return tracts, map_files, map_inputs
