import numpy as np
import healsparse

import lsst.daf.persistence as dafPersist

from .configuration import Configuration
from .patch_mapper import PatchMapper

class TractMapper(object):
    """
    Map a tract.
    """
    def __init__(self, butler, config, outputpath):
        """
        """
        self.butler = butler

        self.config = config

        self.outputpath = outputpath

    def run(self, filter_name, tract):
        """
        """
        skymap = butler.get('deepCoadd_skyMap')

        tract_info = skymap[tract]
        num_patches = tract_info.getNumPatches()

        patch_area = get_sky_polygon_area(tract_info.getPatch_info((0, 0)))
        tract_area = patch_area * num_patches[0] * num_patches[1]
        nside_coverage_tract = 32
        while hp.nside2pixarea(nside_coverage_tract, degrees=True) > tract_area:
            nside_coverage_tract = int(2*nside_coverage_tract)
        nside_coverage_tract = np.clip(nside_coverage_tract / 2, 32, None)

        patch_mapper = PatchMapper(self.butler, self.config, self.outputpath)

        tract_map_list = []
        map_operation_list = []
        for map_type in self.config.map_types.keys():
            op_list = []
            for j, operation in enumerate(self.config.map_types[map_type]):
                op_map_list = []
                tract_map = healsparse.HealSparseMap.make_empty(nside_coverage=nside_coverage_tract,
                                                                nside_sparse=self.config.nside,
                                                                dtype=np.float64)
                op_map_list.append(tract_map)
                op_list.append(op_str_to_code(operation))
            tract_map_list.append(op_map_list)
            map_operation_list.append(op_list)

        for ii in range(num_patches[0]):
            for jj in range(num_patches[1]):
                patch_name = '%d,%d' % (ii, jj)

                patch_input_map, map_values_list = patch_mapper.run(filter_name,
                                                                    tract,
                                                                    patch_name,
                                                                    return_values=True)
                valid_pixels = patch_input_map.valid_pixels

                for i, map_type in enumerate(self.config.map_types.keys()):
                    for j, op in enumerate(map_operation_list[i]):
                        tract_map_list[i][j].update_values_pix(valid_pixels,
                                                               map_values_list[i][j])

        # We are done assembling
        for i, map_type in enumerate(self.config.map_types.keys()):
            for j, op in enumerate(map_operation_list[i]):
                fname = os.path.join(self.outputpath,
                                     self.config.tract_map_filename(filter_name,
                                                                    tract,
                                                                    map_type,
                                                                    op))
                tract_map_list[i][j].write(fname)





