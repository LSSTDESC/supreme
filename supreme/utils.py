import numpy as np
import lsst.sphgeom
import lsst.geom
import lsst.afw.detection as afwDetection
import lsst.meas.algorithms as measAlg


OP_SUM = 1
OP_MEAN = 2
OP_WMEAN = 3
OP_MIN = 4
OP_MAX = 5


def op_code_to_str(op_code):
    if op_code == OP_SUM:
        op_str = 'sum'
    elif op_code == OP_MEAN:
        op_str = 'mean'
    elif op_code == OP_WMEAN:
        op_str = 'wmean'
    elif op_code == OP_MIN:
        op_str = 'min'
    elif op_code == OP_MAX:
        op_str = 'max'

    return op_str


def op_str_to_code(op_str):
    if op_str == 'sum':
        op_code = OP_SUM
    elif op_str == 'mean':
        op_code = OP_MEAN
    elif op_str == 'wmean':
        op_code = OP_WMEAN
    elif op_str == 'min':
        op_code = OP_MIN
    elif op_str == 'max':
        op_code = OP_MAX

    return op_code

def vertices_to_radec(vertices):
    lonlats = [lsst.sphgeom.LonLat(x) for x in vertices]
    radec = np.array([(x.getLon().asDegrees(), x.getLat().asDegrees()) for
                      x in polyLonLats])

    return radec


def approx_patch_polygon_area(patch_info, wcs):
    """
    """
    vertices = patch_info.getInnerSkyPolygon(wcs).getVertices()
    radec = vertices_to_radec(vertices)
    delta_ra = np.max(radec[:, 0]) - np.min(radec[:, 0])
    delta_dec = np.max(radec[:, 1]) - np.min(radec[:, 1])

    # The approximate area is fine for this
    area = delta_ra * delta_dec * np.cos(np.deg2rad(np.mean(radec[:, 1])))

    return area

def radec_to_pixels(wcs, ra, dec):
    pts = [lsst.geom.SpherePoint(r*lsst.geom.degrees, d*lsst.geom.degrees) for
           r, d in zip(ra, dec)]
    pixels = wcs.skyToPixel(pts)

    xy = np.array([(pix.getX(), pix.getY()) for pix in pixels])

    return pixels, xy


def pixels_to_radec(wcs, pixels):
    sph_pts = wcs.pixelToSky(pixels)

    radec = np.array([(sph.getRa().asDegrees(), sph.getDec().asDegrees())
                      for sph in sph_pts])
    return radec


def xy_to_radec(wcs, x, y):
    pixels = [lsst.geom.Point2D(xx, yy) for xx, yy in zip(x, y)]

    return pixels_to_radec(pixels)


def convert_mask_to_bbox_list(mask, plane_name):
    thresh = afwDetection.Threshold(mask.getPlaneBitMask(plane_name),
                                    afwDetection.Threshold.BITMASK)
    fp_list = afwDetection.FootprintSet(mask, thresh).getFootprints()
    defects = measAlg.Defects.fromFootprintList(fp_list)

    return [d.getBBox() for d in defects]


