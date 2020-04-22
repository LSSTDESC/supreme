import os
import numpy as np
import healpy as hp
import healsparse
import esutil

import lsst.geom
import lsst.afw.math as afwMath

from .utils import vertices_to_radec, pixels_to_radec, radec_to_pixels
from .utils import OP_SUM, OP_MEAN, OP_WMEAN, OP_MIN, OP_MAX, OP_OR
from .utils import approx_patch_polygon_area, op_str_to_code
from .utils import convert_mask_to_bbox_list, bbox_to_radec_grid
from .psf import get_approx_psf_size_and_shape


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

        if not self.butler.datasetExists('deepCoadd',
                                         tract=tract,
                                         patch=patch_name,
                                         filter=filter_name):
            if return_values_list:
                return None, None
            else:
                return None

        skymap = self.butler.get('deepCoadd_skyMap')
        tract_info = skymap[tract]

        patch_indices = [int(x) for x in patch_name.split(',')]
        patch_info = tract_info.getPatchInfo(patch_indices)

        # Compute the optimal coverage nside for the size of the patch
        # This does not need to match the tract coverage map!
        patch_area = approx_patch_polygon_area(patch_info, tract_info.getWcs())
        nside_coverage_patch = 32
        while hp.nside2pixarea(nside_coverage_patch, degrees=True) > patch_area:
            nside_coverage_patch = int(2*nside_coverage_patch)
        nside_coverage_patch = int(nside_coverage_patch / 2)

        exposure = self.butler.get('deepCoadd', tract=tract, patch=patch_name, filter=filter_name)
        info = exposure.getInfo()
        inputs = info.getCoaddInputs()

        print('Working on patch %s with %d input ccds' % (patch_name, len(inputs.ccds)))

        # Check if we already have the persisted patch input map
        patch_input_filename = os.path.join(self.outputpath,
                                            self.config.patch_input_filename(filter_name, tract, patch_name))
        if os.path.isfile(patch_input_filename):
            patch_input_map = healsparse.HealSparseMap.read(patch_input_filename)
        else:
            patch_input_map = self.build_patch_input_map(tract_info.getWcs(),
                                                         patch_info,
                                                         inputs.ccds,
                                                         nside_coverage_patch)
            patch_input_map.write(patch_input_filename)

        # Now we build maps
        valid_pixels, ra, dec = patch_input_map.valid_pixels_pos(lonlat=True, return_pixels=True)

        has_psf_quantity = False
        has_metadata_quantity = False
        has_calibrated_quantity = False
        has_coadd_quantity = False
        map_values_list = []
        map_operation_list = []
        for map_type in self.config.map_types.keys():
            if map_type == 'psf_size' or map_type == 'psf_e1' or map_type == 'psf_e2':
                has_psf_quantity = True
            if map_type == 'skylevel' or map_type == 'skysigma' or map_type == 'bgmean':
                has_metadata_quantity = True
            if map_type == 'skylevel' or map_type == 'skysigma' or \
                    map_type == 'bgmean' or map_type == 'background':
                has_calibrated_quantity = True
            if map_type.startswith('coadd'):
                has_coadd_quantity = True

            # Specify maps that are integers not floats
            if map_type in ['nexp', 'coadd_mask']:
                map_dtype = np.int32
            else:
                map_dtype = np.float64

            n_operations = len(self.config.map_types[map_type])
            map_values = np.zeros((valid_pixels.size, n_operations), dtype=map_dtype)
            op_list = []
            for j, operation in enumerate(self.config.map_types[map_type]):
                op_code = op_str_to_code(operation)
                op_list.append(op_code)

                if op_code == OP_MIN or op_code == OP_MAX:
                    # We use fmin and fmax, so nans get overwritten
                    map_values[:, j] = np.nan
                if map_type in ['coadd_image', 'coadd_variance']:
                    if op_code != OP_MEAN:
                        raise RuntimeError("Coadd image, variance must only be MEAN")
                elif map_type in ['coadd_mask']:
                    if op_code != OP_OR:
                        raise RuntimeError("Coadd mask must only be OR")

            map_values_list.append(map_values)
            map_operation_list.append(op_list)

        metadata = patch_input_map.metadata
        weights = np.zeros(valid_pixels.size)
        nexp = np.zeros(valid_pixels.size, dtype=np.int32)
        for bit, ccd in enumerate(inputs.ccds):
            u, = np.where(patch_input_map.check_bits_pix(valid_pixels, [bit]))
            if u.size == 0:
                continue

            # Compute the dataId if we need it (for background)
            dataId = {self.config.detector_id_name: int(ccd['ccd']),
                      self.config.visit_id_name: int(ccd['visit'])}

            # First, the weights and counting, we always need these
            weights[u] += ccd['weight']
            nexp[u] += 1

            if has_psf_quantity:
                ccd_box = lsst.geom.Box2D(ccd.getBBox())
                psf_size, psf_e1, psf_e2 = get_approx_psf_size_and_shape(ccd_box,
                                                                         ccd.getPsf(),
                                                                         ccd.getWcs(),
                                                                         ra[u], dec[u])

            if has_metadata_quantity:
                if ('B%04dSLV' % (bit)) in metadata:
                    skylevel = metadata['B%04dSLV' % (bit)]
                    skysigma = metadata['B%04dSSG' % (bit)]
                    bgmean = metadata['B%04dBGM' % (bit)]
                else:
                    # We must load the metadata
                    # This might be redundant, but useful during testing/development
                    calexp_metadata = self.butler.get('calexp_md', dataId=dataId)
                    skylevel = calexp_metadata['SKYLEVEL']
                    skysigma = calexp_metadata['SKYSIGMA']
                    bgmean = calexp_metadata['BGMEAN']

            if has_calibrated_quantity:
                photoCalib = ccd.getPhotoCalib()
                bf = photoCalib.computeScaledCalibration()
                pixels, xy = radec_to_pixels(ccd.getWcs(), ra[u], dec[u])
                if bf.getBBox() == ccd.getBBox():
                    # This is a variable calibration, track variability
                    calib_scale = photoCalib.getCalibrationMean() * bf.evaluate(xy[:, 0], xy[:, 1])
                else:
                    # This is a spatially constant calibration
                    calib_scale = photoCalib.getCalibrationMean()

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
                elif map_type == 'skylevel':
                    values = skylevel*calib_scale
                elif map_type == 'skysigma':
                    values = skysigma*calib_scale
                elif map_type == 'bgmean':
                    values = bgmean*calib_scale
                elif map_type == 'background':
                    bkg = self.butler.get('calexpBackground', dataId=dataId)
                    bkgImage = bkg.getImage()
                    if self.butler.datasetExists('skyCorr', dataId=dataId):
                        skyCorr = self.butler.get('skyCorr', dataId=dataId)
                        bkgImage += skyCorr.getImage()

                    # Take the background at the given pixel.  Since this
                    # is a smooth map anyway, this should be fine and we
                    # don't need to average over the full coverage
                    values = bkgImage.getArray()[xy[:, 1].astype(np.int32),
                                                 xy[:, 0].astype(np.int32)]
                    values *= calib_scale
                elif map_type.startswith('coadd'):
                    continue
                else:
                    raise ValueError("Illegal map type %s" % (map_type))

                for j, op in enumerate(map_operation_list[i]):
                    if op == OP_SUM:
                        map_values_list[i][u, j] += values
                    elif op == OP_MEAN:
                        map_values_list[i][u, j] += values
                    elif op == OP_WMEAN:
                        map_values_list[i][u, j] += ccd['weight']*values
                    elif op == OP_MIN:
                        map_values_list[i][u, j] = np.fmin(map_values_list[i][u, j], values)
                    elif op == OP_MAX:
                        map_values_list[i][u, j] = np.fmax(map_values_list[i][u, j], values)
                    elif op == OP_OR:
                        # This is not actually supported by anything yet
                        map_values_list[i][u, j] = np.bitwise_or(map_values_list[i][u, j],
                                                                 values.astype(np.int64))

        if has_coadd_quantity:
            # Convert the coadd pixel grid to ra/dec
            coadd_xy, coadd_radec = bbox_to_radec_grid(exposure.getWcs(),
                                                       patch_info.getInnerBBox())
            coadd_origin = exposure.getBBox().getBegin()
            coadd_xy[:, 0] -= coadd_origin.getX()
            coadd_xy[:, 1] -= coadd_origin.getY()

            # Convert radec to healpix pixels
            hpix = hp.ang2pix(self.config.nside, coadd_radec[:, 0], coadd_radec[:, 1],
                              lonlat=True, nest=True)
            minpix = np.min(hpix)
            maxpix = np.max(hpix)

            # Need to map the coadd pixels to the valid pixels
            aa, bb = esutil.numpy_util.match(valid_pixels - minpix, np.arange(maxpix - minpix + 1))

            npix_arr = np.zeros(maxpix - minpix + 1, dtype=np.int32)

            np.add.at(npix_arr,
                      hpix - minpix,
                      1)
            coadd_pix_use, = np.where(npix_arr > 0)

            for i, map_type in enumerate(self.config.map_types.keys()):
                if map_type == 'coadd_image':
                    coadd_array = exposure.getImage().getArray()
                elif map_type == 'coadd_variance':
                    coadd_array = exposure.getVariance().getArray()
                elif map_type == 'coadd_mask':
                    coadd_array = exposure.getMask().getArray()
                else:
                    # Not a coadd map, skip this.
                    continue

                for j, op in enumerate(map_operation_list[i]):
                    if op == OP_SUM or op == OP_MEAN:
                        values_arr = np.zeros_like(npix_arr, dtype=np.float64)
                        np.add.at(values_arr,
                                  hpix - minpix,
                                  coadd_array[coadd_xy[:, 1], coadd_xy[:, 0]])
                        if op == OP_MEAN:
                            values_arr[coadd_pix_use] /= npix_arr[coadd_pix_use]
                    elif op == OP_MIN:
                        values_arr = np.zeros_like(npix_arr, dtype=np.float64) + np.nan
                        np.fmin.at(values_arr,
                                   hpix - minpix,
                                   coadd_array[coadd_xy[:, 1], coadd_xy[:, 0]])
                    elif op == OP_MAX:
                        values_arr = np.zeros_like(npix_arr, dtype=np.float64) + np.nan
                        np.fmax.at(values_arr,
                                   hpix - minpix,
                                   coadd_array[coadd_xy[:, 1], coadd_xy[:, 0]])
                    elif op == OP_OR:
                        values_arr = np.zeros_like(npix_arr, dtype=np.int32)
                        np.bitwise_or.at(values_arr,
                                         hpix - minpix,
                                         coadd_array[coadd_xy[:, 1], coadd_xy[:, 0]])

                # And convert these back to the map_values array
                map_values_list[i][aa, j] = values_arr[bb]

        # And we've done all the accumulations, finish the mean/wmean
        # And turn these into maps
        for i, map_type in enumerate(self.config.map_types.keys()):
            for j, op in enumerate(map_operation_list[i]):
                if not map_type.startswith('coadd'):
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
                                                                   dtype=map_values_list[i][:, j].dtype)
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

            dataId = {self.config.detector_id_name: int(ccd['ccd']),
                      self.config.visit_id_name: int(ccd['visit'])}

            if self.config.use_calexp_mask:
                calexp = self.butler.get('calexp', dataId=dataId)
                calexp_metadata = calexp.getMetadata()

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
                                                                   dtype=np.uint8)
                    healsparse.realize_geom(mask_poly_list, mask_map)
                    poly_map.apply_mask(mask_map)

                if 'SKYLEVEL' not in calexp_metadata:
                    # We must recompute skylevel, skysigma

                    # We need to re-add in the background
                    bkg = self.butler.get('calexpBackground', dataId=dataId)

                    statsControl = afwMath.StatisticsControl(3.0, 3)
                    maskVal = calexp.getMaskedImage().getMask().getPlaneBitMask(["BAD",
                                                                                 "SAT",
                                                                                 "DETECTED"])
                    statsControl.setAndMask(maskVal)
                    maskedImage = calexp.getMaskedImage()
                    maskedImage += bkg.getImage()
                    stats = afwMath.makeStatistics(maskedImage, afwMath.MEDIAN | afwMath.STDEVCLIP,
                                                   statsControl)
                    skylevel = stats.getValue(afwMath.MEDIAN)
                    skysigma = stats.getValue(afwMath.STDEVCLIP)
                    del maskedImage
                else:
                    skylevel = calexp_metadata['SKYLEVEL']
                    skysigma = calexp_metadata['SKYSIGMA']
            else:
                calexp_metadata = self.butler.get('calexp_md', dataId=dataId)
                if 'SKYLEVEL' not in calexp_metadata:
                    # We want to log this
                    skylevel = 0.0
                    skysigma = 0.0
                else:
                    skylevel = calexp_metadata['SKYLEVEL']
                    skysigma = calexp_metadata['SKYSIGMA']

            metadata['B%04dSLV' % (bit)] = skylevel
            metadata['B%04dSSG' % (bit)] = skysigma
            metadata['B%04dBGM' % (bit)] = calexp_metadata['BGMEAN']
            metadata['B%04dBGV' % (bit)] = calexp_metadata['BGVAR']

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
