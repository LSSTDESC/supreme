import os
import numpy as np
import healsparse

import lsst.daf.persistence as dafPersist

from .configuration import Configuration
from .utils import vertices_to_radec, radec_to_pixels, pixels_to_radec, xy_to_radec
from .utils import OP_SUM, OP_MEAN, OP_WMEAN, OP_MIN, OP_MAX

class PatchMapper(object):
    """
    Map a patch.
    """
    def __init__(self, butler, config, outputpath):
        """
        """
        self.butler = butler

        self.config = config

        self.outputpath = outputpath

    def run(self, filter_name, tract, patch_name, return_values_list=True):
        """
        """
        print('Working on patch %s' % (patch_name))

        if not self.butler.datasetExists('deepCoadd',
                                         tract=tract,
                                         patch=patch_name,
                                         filter=filter_name):
            return

        skymap = butler.get('deepCoadd_skyMap')
        tract_info = skymap[tract]

        patch_indices = [int(x) for x in patch_name.split(',')]
        patch_info = tractInfo.getPatchInfo(patch_indices)

        # Compute the optimal coverage nside for the size of the patch
        # This does not need to match the tract coverage map!
        patch_area = get_sky_polygon_area(patch_info, tract_info.getWcs())
        nside_coverage_patch = 32
        while hp.nside2pixarea(nside_coverage_patch, degrees=True) > patch_area:
            nside_coverage_patch = int(2*nside_coverage_patch)
        nside_coverage_patch = int(nside_coverage_patch / 2)

        exposure = butler.get('deepCoadd', tract=tract, patch=patch_name, filter=filter_name)
        info = exposure.getInfo()
        inputs = info.getCoaddInputs()

        # Check if we already have the persisted patch input map
        patch_input_filename = os.path.join(self.outputpath,
                                            self.config.patch_input_filename(filter_name, tract, patch_name))
        if os.path.isfile(patch_input_filename):
            patch_input_map = healsparse.HealSparseMap.read(patch_input_filename)
        else:
            patch_input_map = self.build_patch_input_map(tract_info.getWcs(), patch_info, ccds, nside_coverage_patch)
            patch_input_map.write(patch_input_filename)

        # Now we build maps
        valid_pixels, ra, dec = patch_input_map.valid_pixels_pos(lonlat=True, return_pixels=True)

        has_psf_quantity = False
        save_nexp = False
        map_value_list = []
        map_operation_list = []
        for map_type in self.config.map_types.keys():
            if map_type == 'psf_size' or map_type == 'psf_e1' or map_type == 'psf_e2':
                has_psf_quantity = True
            elif map_type == 'nexp':
                save_nexp = True

            n_operations = len(self.config.map_types[map_type])
            map_values = np.zeros((valid_pixels.size, n_operations))
            op_list = []
            for j, operation in enumerate(self.config.map_types[map_type]):
                op_code = op_str_to_code(operation)
                op_list.append(op_code)

            map_value_list.append(map_values)
            map_operation_list.append(op_list)

        weights = np.zeros(valid_pixels.size)
        nexp = np.zeros(valid_pixels.size, dtype=np.int32)
        for bit, ccd in enumerate(inputs.ccds):
            u, = np.where(patch_input_map.check_bits_pix(valid_pixels, [bit]))
            if u.size == 0:
                continue

            # First, the weights and counting, we always need these
            weights[u] += ccd['weight']
            nexp[u] += 1

            if has_psf_quantities:
                ccd_box = lsst.geom.Box2D(ccd.getBBox())
                psf_size, psf_e1, psf_e2 = get_approx_psf_size_and_shape(ccd_box,
                                                                         ccd.getPsf(),
                                                                         ccd.getWcs(),
                                                                         ra[u], dec[u])

            for i, map_type in enumerate(self.config.map_types.keys()):
                if map_type == 'psf_size':
                    values = psf_size
                elif map_type == 'psf_e1':
                    values = psf_e1
                elif map_type == 'psf_e2':
                    values = psf_e2
                elif map_type == 'exptime':
                    values = np.zeros(u.size) + ccd.getVisitInfo().getExposureTime()
                elif map_type == 'airmass':
                    values = np.zeros(u.size) + ccd.getVisitInfo().getBoresightAirmass()
                elif map_type == 'nexp':
                    values = np.ones(u.size, dtype=np.int32)
                else:
                    raise ValueError("Illegal map type %s" % (map_type))

                for j, op in enumerate(map_operation_list[i]):
                    if op == OP_SUM:
                        map_values_list[i][u, j] += values
                    elif op == OP_MEAN:
                        map_values_list[i][u, j] += values
                    elif op == OP_WMEAN:
                        map_values_list[i][u, j] += weights[u]*values
                    elif op == OP_MIN:
                        map_values_list[i][u, j] = np.fmin(map_values_list[i][u, j], values)
                    elif op == OP_MAX:
                        map_values_list[i][u, j] = np.fmax(map_values_list[i][u, j], values)

        # And we've done all the accumulations, finish the mean/wmean
        # And turn these into maps
        returned_maps = []
        for i, map_type in enumerate(self.config.map_types.keys()):
            for j, op in enumerate(map_operation_list[i]):
                if op == OP_MEAN:
                    map_values_list[i][:, j] /= nexp
                elif op == OP_WMEAN:
                    map_values_list[i][:, j] /= weights

                if not return_values_list:
                    # Here we want to save the patch map
                    fname = os.path.join(self.outputpath,
                                         self.config.patch_map_filename(filter_name,
                                                                        tract,
                                                                        patch_name,
                                                                        map_type,
                                                                        op))
                    temp_map = healsparse.HealSparseMap.make_empty(nside_coverage=nside_coverage_patch,
                                                                   nside_sparse=self.config.nside,
                                                                   dtype=np.float64)
                    temp_map.update_values_pix(valid_pixels, map_values_list[i][:, j])
                    temp_map.write(fname)

        if return_values_list:
            return patch_input_map, map_values_list

    def build_patch_input_map(self, tract_wcs, patch_info, ccds, nside_coverage_patch):
        """
        """
        patch_input_map = healsparse.HealSparseMap.make_empty(nside_coverage=nside_coverage_patch,
                                                              nside_sparse=self.config.nside,
                                                              dtype=healsparse.WIDE_MASK,
                                                              wide_mask_maxbits=len(ccds))
        metadata = {}
        poly_map_list = []
        for bit, ccd in enumerate(ccds):
            metadata['B%04dCCD' % (bit)] = ccd['ccd']
            metadata['B%04dVIS' % (bit)] = ccd['visit']
            metadata['B%04dWT' % (bit)] = ccd['weight']

            wcs = ccd.getWcs()
            ccd_poly = ccd.getValidPolygon().convexHull()
            ccd_poly_radec = pixels_to_radec(wcs, ccd_poly.getVertices())

            # Use polygons for all of these
            poly = healsparse.Polygon(ra=ccd_poly_radec[: -1, 0],
                                      dec=ccd_poly_radec[: -1, 1],
                                      value=[bit])
            poly_map = poly.get_map_like(patch_input_map)

            if self.config.use_calexp_mask:
                calexp = self.butler.get('calexp', ccd=ccd['ccd'], visit=ccd['visit'])
                mask = calexp.getMask()

                mask_poly_list = []
                for plane in self.config.bad_mask_planes:
                    bboxes = convert_mask_to_bbox_list(mask, plane)
                    for bbox in bboxes:
                        b_radec = pixels_to_radec(wcs, lsst.geom.Box2D(bbox).getCorners())
                        mask_poly = healsparse.Polygon(ra=b_radec[:, 0],
                                                       dec=b_radec[:, 1],
                                                       value=1)
                        mask_poly_list.append(mask_poly)
                    mask_map = healsparse.HealSparseMap.make_empty(nside_coverage=nside_coverage_patch,
                                                                   nside_sparse=self.config.nside,
                                                                   dtype=uint8)
                    healsparse.realize_geom(mask_poly_list, mask_map)

                    # valid_pixels = poly_map.valid_pixels
                    # bad, = np.where(mask_map.get_values_pix(valid_pixels) == 1)
                    # poly_map.clear_bits_pix(valid_pixels[bad], [bit])
                    poly_map.apply_mask(mask_map)

            # Now we have the full masked ccd map, append to list
            poly_map_list.append(poly_map)

        # Now cut down to the inner tract polygon
        poly_vertices = patch_info.getInnerSkyPolygon(tract_wcs).getVertices()
        patch_radec = vertices_to_radec(poly_vertices)
        patch_poly = healsparse.Polygon(ra=patch_radec[:, 0], dec=patch_radec[:, 1],
                                        value=np.arange(patch_input_map.wide_mask_maxbits))

        # Realize the patch polygon
        patch_poly_map = patch_poly.get_map_like(patch_input_map)
        # Combine all masked ccd polygons
        patch_input_map = healsparse.or_union(poly_map_list)
        # And limit to the patch polygon
        patch_input_map = healsparse.and_intersection([patch_input_map, patch_poly_map])

        # And set the metadata
        patch_input_map.metadata = metadata

        return patch_input_map
