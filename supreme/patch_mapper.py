import os
import multiprocessing
import numpy as np
import healpy as hp
import healsparse
import esutil
import warnings

import lsst.geom
import lsst.afw.geom as afwGeom
import lsst.afw.math as afwMath
import lsst.daf.persistence as dafPersist

import astropy.units as units
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz

from .utils import vertices_to_radec, pixels_to_radec, radec_to_xy
from .utils import OP_NONE, OP_SUM, OP_MEAN, OP_WMEAN, OP_MIN, OP_MAX, OP_OR
from .utils import approx_patch_polygon_area, op_str_to_code
from .utils import convert_mask_to_bbox_list, bbox_to_radec_grid
from .psf import get_approx_psf_size_and_shape

global _butler_dict


_butler_dict = {}


def pool_initializer(butler, is_single):
    proc = multiprocessing.Process()
    this_core = proc._identity[0]

    if is_single:
        _butler_dict[this_core] = butler
    else:
        repo = butler._repos._inputs[0].repoArgs.root
        _butler_dict[this_core] = dafPersist.Butler(repo)


class PatchMapper(object):
    """
    Map a single patch.  Should usually be invoked via MultiMapper.
    """
    def __init__(self, outputpath, config):
        """
        Instantiate a PatchMapper

        Parameters
        ----------
        outputpath : `str`
           Base path for the output files.
        config : `Configuration`
           supreme configuration object

        Returns
        -------
        patchMapper : `supreme.PatchMapper`
        """
        self.outputpath = outputpath
        self.config = config

        if not os.path.isdir(outputpath):
            raise RuntimeError("Outputpath %s does not exist." % (outputpath))

    def __call__(self, tract, filter_name, patch_name, tract_mode=True,
                 map_run=None, clobber=False):
        """
        Compute the map for a single patch.

        Parameters
        ----------
        tract : `int`
           Tract number
        filter_name : `str`
           Filter name (dataId format)
        patch_name : `str`
           Patch name
        tract_mode : `bool`, optional
           If being run in tract_mode, then the values will be returned as a list
           of arrays.  If False, the patch map will be written to disk.
        map_run : `dict` or None, optional
           Dict matching config.map_types, with a True/False flag if a map
           is to be run (for use when a map already exists).  Default is None,
           which means run all.
        clobber : `bool`, optional
           Clobber any existing files.

        Returns
        -------
        patch_map : `healsparse.HealSparseMap`
           Sparse "input" map encoding the coverage of the patch.
        values : `list` of `numpy.ndarray`, optional
           List of value arrays for each map/operation type specified in the config.
           Returned if return_values_list is True.
        """
        # Copy in global butler for multiprocessing
        this_core = multiprocessing.current_process()._identity[0]

        butler = _butler_dict[this_core]

        if not butler.datasetExists('deepCoadd',
                                    tract=tract,
                                    patch=patch_name,
                                    filter=filter_name):
            if tract_mode:
                return None, None, False
            else:
                return None, False

        # Create the path for tract/patches (if nessary)
        os.makedirs(os.path.join(self.outputpath,
                                 self.config.patch_relpath(tract)),
                    exist_ok=True)

        skymap = butler.get('deepCoadd_skyMap')
        tract_info = skymap[tract]

        patch_indices = [int(x) for x in patch_name.split(',')]
        patch_info = tract_info.getPatchInfo(patch_indices)

        nside_coverage_patch = self._compute_nside_coverage_patch(patch_info, tract_info)
        try:
            exposure = butler.get('deepCoadd', tract=tract, patch=patch_name, filter=filter_name)
        except RuntimeError:
            print("Error reading deepCoadd for %d / %s / %s" % (tract, filter_name, patch_name))
            if tract_mode:
                return None, None, True
            else:
                return None, True

        info = exposure.getInfo()
        inputs = info.getCoaddInputs()

        print('Working on patch %s with %d input ccds' % (patch_name, len(inputs.ccds)))

        # Check if we already have the persisted patch input map
        patch_input_filename = os.path.join(self.outputpath,
                                            self.config.patch_relpath(tract),
                                            self.config.patch_input_filename(filter_name, tract, patch_name))
        if os.path.isfile(patch_input_filename) and not clobber:
            patch_input_map = healsparse.HealSparseMap.read(patch_input_filename)
        else:
            patch_input_map = self.build_patch_input_map(butler,
                                                         tract_info.getWcs(),
                                                         patch_info,
                                                         inputs.ccds,
                                                         nside_coverage_patch)
            patch_input_map.write(patch_input_filename, clobber=clobber)

        # Now we build maps
        valid_pixels, vpix_ra, vpix_dec = patch_input_map.valid_pixels_pos(lonlat=True,
                                                                           return_pixels=True)

        has_psf_quantity = False
        has_metadata_quantity = False
        has_calibrated_quantity = False
        has_coadd_quantity = False
        has_parallactic_quantity = False
        has_zenith_quantity = False
        map_values_list = []
        self.map_operation_list = []
        for map_type in self.config.map_types.keys():
            # Check if we need to run _any_ operations
            if map_run is None or np.any(map_run[map_type]):
                if map_type == 'psf_size' or map_type == 'psf_e1' or map_type == 'psf_e2':
                    has_psf_quantity = True
                if map_type == 'skylevel' or map_type == 'skysigma' or map_type == 'bgmean':
                    has_metadata_quantity = True
                if map_type == 'skylevel' or map_type == 'skysigma' or \
                        map_type == 'bgmean' or map_type == 'background':
                    has_calibrated_quantity = True
                if map_type == 'airmass':
                    has_zenith_quantity = True
                if map_type.startswith('dcr') or map_type == 'parallactic':
                    has_parallactic_quantity = True
                    has_zenith_quantity = True
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
                if map_run is not None and not map_run[map_type][j]:
                    op_code = OP_NONE
                op_list.append(op_code)

                if op_code == OP_MIN or op_code == OP_MAX:
                    # We use fmin and fmax, so nans get overwritten
                    map_values[:, j] = np.nan
                if map_type in ['coadd_image', 'coadd_variance']:
                    if op_code != OP_MEAN and op_code != OP_NONE:
                        raise RuntimeError("Coadd image, variance must only be MEAN")
                elif map_type in ['coadd_mask']:
                    if op_code != OP_OR and op_code != OP_NONE:
                        raise RuntimeError("Coadd mask must only be OR")

            map_values_list.append(map_values)
            self.map_operation_list.append(op_list)

        metadata = patch_input_map.metadata
        weights = np.zeros(valid_pixels.size)
        nexp = np.zeros(valid_pixels.size, dtype=np.int32)
        loc = None
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

            if has_zenith_quantity:
                if loc is None:
                    # Get the observatory location once
                    obs = ccd.getVisitInfo().getObservatory()
                    loc = EarthLocation(lat=obs.getLatitude().asDegrees()*units.deg,
                                        lon=obs.getLongitude().asDegrees()*units.deg,
                                        height=obs.getElevation()*units.m)

                zenith = self._compute_zenith_angle(loc, ccd.getVisitInfo().getDate().get(),
                                                    np.median(vpix_ra[u]),
                                                    np.median(vpix_dec[u]))

            if has_psf_quantity:
                ccd_box = lsst.geom.Box2D(ccd.getBBox())
                psf_size, psf_e1, psf_e2 = get_approx_psf_size_and_shape(ccd_box,
                                                                         ccd.getPsf(),
                                                                         ccd.getWcs(),
                                                                         vpix_ra[u],
                                                                         vpix_dec[u])

            if has_metadata_quantity:
                if ('B%04dSLV' % (bit)) in metadata:
                    skylevel = metadata['B%04dSLV' % (bit)]
                    skysigma = metadata['B%04dSSG' % (bit)]
                    bgmean = metadata['B%04dBGM' % (bit)]
                else:
                    # We must load the metadata
                    # This might be redundant, but useful during testing/development
                    calexp_metadata = butler.get('calexp_md', dataId=dataId)
                    skylevel = calexp_metadata['SKYLEVEL']
                    skysigma = calexp_metadata['SKYSIGMA']
                    bgmean = calexp_metadata['BGMEAN']

            if has_parallactic_quantity:
                vi = ccd.getVisitInfo()
                par_angle = vi.getBoresightParAngle().asRadians()

            if has_calibrated_quantity:
                xy = radec_to_xy(ccd.getWcs(), vpix_ra[u], vpix_dec[u])
                calib_scale = self._compute_calib_scale(ccd, xy)

            for i, map_type in enumerate(self.config.map_types.keys()):
                if map_run is not None and not np.any(map_run[map_type]):
                    values = 0
                elif map_type == 'psf_size':
                    values = psf_size
                elif map_type == 'psf_e1':
                    values = psf_e1
                elif map_type == 'psf_e2':
                    values = psf_e2
                elif map_type == 'exptime':
                    values = np.zeros(u.size) + ccd.getVisitInfo().getExposureTime()
                elif map_type == 'airmass':
                    # Use the median position on the CCD for the airmass
                    values = self._compute_airmass(zenith)
                elif map_type == 'boresight_dist':
                    # Distance from the boresight in radians
                    bore = ccd.getVisitInfo().getBoresightRaDec()
                    bore_ra_rad, bore_dec_rad = bore.getRa().asRadians(), bore.getDec().asRadians()
                    values = esutil.coords.sphdist(bore_ra_rad, bore_dec_rad,
                                                   np.deg2rad(vpix_ra[u]),
                                                   np.deg2rad(vpix_dec[u]),
                                                   units=['rad', 'rad'])
                elif map_type == 'dcr_dra':
                    values = np.zeros(u.size) + np.tan(zenith)*np.sin(par_angle)
                elif map_type == 'dcr_ddec':
                    values = np.zeros(u.size) + np.tan(zenith)*np.cos(par_angle)
                elif map_type == 'dcr_e1':
                    values = np.zeros(u.size) + (np.tan(zenith)**2.)*np.cos(2*par_angle)
                elif map_type == 'dcr_e2':
                    values = np.zeros(u.size) + (np.tan(zenith)**2.)*np.sin(2*par_angle)
                elif map_type == 'parallactic':
                    values = np.zeros(u.size) + par_angle
                elif map_type == 'nexp':
                    values = np.ones(u.size, dtype=np.int32)
                elif map_type == 'skylevel':
                    values = skylevel*calib_scale
                elif map_type == 'skysigma':
                    values = skysigma*calib_scale
                elif map_type == 'bgmean':
                    values = bgmean*calib_scale
                elif map_type == 'background':
                    values = self._compute_calexp_background(butler, dataId, xy)
                    values *= calib_scale
                elif map_type.startswith('coadd'):
                    continue
                else:
                    raise ValueError("Illegal map type %s" % (map_type))

                for j, op in enumerate(self.map_operation_list[i]):
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
            self._update_coadd_map_values(exposure, patch_info, valid_pixels, map_values_list)

        # And we've done all the accumulations, finish the mean/wmean
        # And turn these into maps
        for i, map_type in enumerate(self.config.map_types.keys()):
            for j, op in enumerate(self.map_operation_list[i]):
                if not map_type.startswith('coadd'):
                    if op == OP_MEAN:
                        map_values_list[i][:, j] /= nexp
                    elif op == OP_WMEAN:
                        map_values_list[i][:, j] /= weights

                if not tract_mode:
                    # Here we want to save the patch map
                    fname = os.path.join(self.outputpath,
                                         self.config.patch_relpath(tract),
                                         self.config.patch_map_filename(filter_name,
                                                                        tract,
                                                                        patch_name,
                                                                        map_type,
                                                                        op))
                    temp_map = healsparse.HealSparseMap.make_empty(nside_coverage=nside_coverage_patch,
                                                                   nside_sparse=self.config.nside,
                                                                   dtype=map_values_list[i][:, j].dtype)
                    temp_map.update_values_pix(valid_pixels, map_values_list[i][:, j])
                    if not os.path.isfile(fname) or clobber:
                        # Only write if clobber is True or if file does not exist
                        temp_map.write(fname, clobber=clobber)

        if tract_mode:
            return patch_input_map, map_values_list, False
        else:
            return None, False

    def build_patch_input_map(self, butler, tract_wcs, patch_info, ccds, nside_coverage_patch):
        """
        Build the patch input map.

        Parameters
        ----------
        butler : `lsst.daf.persistence.Butler`
           gen2 butler
        tract_wcs : `lsst.afw.geom.SkyWcs`
           WCS object for the tract
        patch_info : `lsst.skymap.PatchInfo`
           Patch info object
        ccds : `lsst.afw.table.ExposureCatalog`
           Catalog of ccd information
        nside_coverage_patch : `int`
           Healpix nside for the coverage map

        Returns
        -------
        patch_input_map : `healsparse.HealSparseMap`
           Healsparse map encoding input ccd information.
        """
        patch_input_map = healsparse.HealSparseMap.make_empty(nside_coverage=nside_coverage_patch,
                                                              nside_sparse=self.config.nside,
                                                              dtype=healsparse.WIDE_MASK,
                                                              wide_mask_maxbits=len(ccds))
        metadata = {}
        for bit, ccd in enumerate(ccds):
            metadata['B%04dCCD' % (bit)] = ccd['ccd']
            metadata['B%04dVIS' % (bit)] = ccd['visit']
            metadata['B%04dWT' % (bit)] = ccd['weight']

            wcs = ccd.getWcs()
            ccd_poly = ccd.getValidPolygon()
            if ccd_poly is None:
                ccd_poly = afwGeom.Polygon(lsst.geom.Box2D(ccd.getBBox()))
            ccd_poly_radec = pixels_to_radec(wcs, ccd_poly.convexHull().getVertices())

            # Use polygons for all of these
            poly = healsparse.Polygon(ra=ccd_poly_radec[: -1, 0],
                                      dec=ccd_poly_radec[: -1, 1],
                                      value=[bit])
            poly_map = poly.get_map_like(patch_input_map)

            dataId = {self.config.detector_id_name: int(ccd['ccd']),
                      self.config.visit_id_name: int(ccd['visit'])}

            if self.config.use_calexp_mask:
                calexp = butler.get('calexp', dataId=dataId)
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
                    skylevel, skysigma = self._compute_skylevel(butler, dataId, calexp)
                else:
                    skylevel = calexp_metadata['SKYLEVEL']
                    skysigma = calexp_metadata['SKYSIGMA']
            else:
                calexp_metadata = butler.get('calexp_md', dataId=dataId)
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

            # Now we have the full masked ccd map, set the appropriate bit
            patch_input_map.set_bits_pix(poly_map.valid_pixels, [bit])

        # Now cut down to the inner tract polygon
        poly_vertices = patch_info.getInnerSkyPolygon(tract_wcs).getVertices()
        patch_radec = vertices_to_radec(poly_vertices)
        patch_poly = healsparse.Polygon(ra=patch_radec[:, 0], dec=patch_radec[:, 1],
                                        value=np.arange(patch_input_map.wide_mask_maxbits))

        # Realize the patch polygon
        patch_poly_map = patch_poly.get_map_like(patch_input_map)
        patch_input_map = healsparse.and_intersection([patch_input_map, patch_poly_map])

        # And set the metadata
        patch_input_map.metadata = metadata

        return patch_input_map

    def _compute_skylevel(self, butler, dataId, calexp):
        """
        Compute the skylevel and skysigma for a calexp.

        Parameters
        ----------
        butler : `lsst.daf.persistence.Butler`
           gen2 butler
        dataId : `dict`
           dataId dictionary for the individual calexp.
        calexp : `lsst.afw.image.ExposureF`
           Exposure to measure sky level.

        Returns
        -------
        skylevel : `float`
           Level of the sky associated with calexp
        skysigma : `float`
           Stddev of sky associated with calexp
        """
        bkg = butler.get('calexpBackground', dataId=dataId)

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

        return skylevel, skysigma

    def _compute_calib_scale(self, ccd, xy):
        """
        Compute calibration scaling value

        Parameters
        ----------
        ccd : `lsst.afw.table.ExposureRecord`
           Exposure metadata for given ccd
        xy : `numpy.ndarray`
           Nx2 array of x/y positions to compute calibration scale.

        Returns
        -------
        calib_scale : `numpy.ndarray`
           Length N array of calibration scale values
        """
        photoCalib = ccd.getPhotoCalib()
        bf = photoCalib.computeScaledCalibration()
        if bf.getBBox() == ccd.getBBox():
            # This is a variable calibration, track variability
            calib_scale = photoCalib.getCalibrationMean() * bf.evaluate(xy[:, 0], xy[:, 1])
        else:
            # This is a spatially constant calibration
            calib_scale = photoCalib.getCalibrationMean()

        return calib_scale

    def _compute_calexp_background(self, butler, dataId, xy):
        """
        Compute background value for a calexp at a list of positions.  Uses
        calexpBackground and skyCorr (if available).

        Parameters
        ----------
        butler : `lsst.daf.persistence.Butler`
           gen2 butler
        dataId : `dict`
           dataId dictionary for the individual calexp.
        xy : `numpy.ndarray`
           Nx2 array of x/y positions to compute calibration scale.

        Returns
        -------
        values : `numpy.ndarray`
           Length N array of background values.
        """
        bkg = butler.get('calexpBackground', dataId=dataId)
        bkgImage = bkg.getImage()
        if butler.datasetExists('skyCorr', dataId=dataId):
            skyCorr = butler.get('skyCorr', dataId=dataId)
            bkgImage += skyCorr.getImage()

        # Take the background at the given pixel.  Since this
        # is a smooth map anyway, this should be fine and we
        # don't need to average over the full coverage
        values = bkgImage.getArray()[xy[:, 1].astype(np.int32),
                                     xy[:, 0].astype(np.int32)]

        return values

    def _compute_zenith_angle(self, loc, mjd, ra, dec):
        """
        Compute the zenith angle for one or many ra/dec.

        Parameters
        ----------
        loc : `astropy.coordinates.earth.EarthLocation`
        mjd : `float`
        ra : `np.ndarray`
           Right ascension
        dec : `np.ndarray`
           Declination

        Returns
        -------
        zenith : `np.ndarray`
           Zenith angle(s) in radians
        """
        t = Time(mjd, format='mjd')
        c = SkyCoord(ra, dec, unit='deg')
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            c_altaz = c.transform_to(AltAz(obstime=t, location=loc))
        return np.pi/2. - c_altaz.alt.rad

    def _compute_airmass(self, zenith):
        """
        Compute the airmass for a list of zenith angles.
        Computed using simple expansion formula.

        Parameters
        ----------
        zenith : `np.ndarray`
           Zenith angle(s), radians

        Returns
        -------
        airmass : `np.ndarray`
        """
        secz = 1./np.cos(zenith)
        airmass = (secz -
                   0.0018167*(secz - 1.0) -
                   0.002875*(secz - 1.0)**2.0 -
                   0.0008083*(secz - 1.0)**3.0)
        return airmass

    def _update_coadd_map_values(self, exposure, patch_info, valid_pixels, map_values_list):
        """
        Modify the map values in-place for coadd quantities

        Parameters
        ----------
        exposure : `lsst.afw.image.ExposureF`
           Coadd exposure object
        patch_info : `lsst.afw.table.ExposureRecord`
           Patch information
        valid_pixels : `numpy.ndarray`
           Indices of valid pixels to compute values
        map_values_list : `list`
           List encoding all the map values to accumulate.
        """
        if valid_pixels.size == 0:
            # There is nothing to do here.
            return

        # Convert the coadd pixel grid to ra/dec
        coadd_xy, coadd_radec = bbox_to_radec_grid(exposure.getWcs(),
                                                   patch_info.getInnerBBox())
        coadd_origin = exposure.getBBox().getBegin()
        coadd_xy[:, 0] -= coadd_origin.getX()
        coadd_xy[:, 1] -= coadd_origin.getY()

        # Apply the mask, to remove pixels with no coadd data
        mask = exposure.getMask()
        mask_array = mask.getArray()
        good, = np.where((mask_array[coadd_xy[:, 1], coadd_xy[:, 0]] &
                          2**mask.getMaskPlaneDict()['NO_DATA']) == 0)

        if good.size == 0:
            # There are no good pixels with data, so just kick out
            return

        coadd_xy = coadd_xy[good, :]
        coadd_radec = coadd_radec[good, :]

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

            for j, op in enumerate(self.map_operation_list[i]):
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

    def _compute_nside_coverage_patch(self, patch_info, tract_info):
        """
        Compute the optimal coverage nside for a patch.

        Parameters
        ----------
        patch_info : `lsst.skymap.PatchInfo`
           Information on the patch
        tract_info : `lsst.skymap.ExplicitTractInfo`
           Information on the tract

        Returns
        -------
        nside_coverage_patch : `int`
           Optimal coverage nside
        """
        # Compute the optimal coverage nside for the size of the patch
        # This does not need to match the tract coverage map!
        patch_area = approx_patch_polygon_area(patch_info, tract_info.getWcs())
        nside_coverage_patch = 32
        while hp.nside2pixarea(nside_coverage_patch, degrees=True) > patch_area:
            nside_coverage_patch = int(2*nside_coverage_patch)
        nside_coverage_patch = int(nside_coverage_patch / 2)

        return nside_coverage_patch
