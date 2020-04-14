
import numpy as np
import lsst.geom
from lsst.afw.math import ChebyshevBoundedField, ChebyshevBoundedFieldControl


def get_approx_psf_size_and_shape(bbox, psf, wcs, ra, dec, nx=20, ny=20, orderx=2, ordery=2):
    pts = [lsst.geom.SpherePoint(r*lsst.geom.degrees, d*lsst.geom.degrees) for
           r, d in zip(ra, dec)]
    pixels = wcs.skyToPixel(pts)

    ctrl = ChebyshevBoundedFieldControl()
    ctrl.orderX = orderx
    ctrl.orderY = ordery
    ctrl.triangular = False

    xSteps = np.linspace(bbox.getMinX(), bbox.getMaxX(), nx)
    ySteps = np.linspace(bbox.getMinY(), bbox.getMaxY(), ny)
    x = np.tile(xSteps, nx)
    y = np.repeat(ySteps, ny)

    psf_size = np.zeros(x.size)
    psf_e1 = np.zeros(x.size)
    psf_e2 = np.zeros(x.size)

    for i in range(x.size):
        shape = psf.computeShape(lsst.geom.Point2D(x[i], y[i]))
        psf_size[i] = shape.getDeterminantRadius()
        ixx = shape.getIxx()
        iyy = shape.getIyy()
        ixy = shape.getIxy()

        psf_e1[i] = (ixx - iyy)/(ixx + iyy + 2.*psf_size[i]**2.)
        psf_e2[i] = (2.*ixy)/(ixx + iyy + 2.*psf_size[i]**2.)

    pixel_x = np.array([pix.getX() for pix in pixels])
    pixel_y = np.array([pix.getY() for pix in pixels])

    cheb_size = ChebyshevBoundedField.fit(lsst.geom.Box2I(bbox), x, y, psf_size, ctrl)
    psf_size_pts = chebSize.evaluate(pixel_x, pixel_y)
    cheb_e1 = ChebyshevBoundedField.fit(lsst.geom.Box2I(bbox), x, y, psf_e1, ctrl)
    psf_e1_pts = cheb_e1.evaluate(pixel_x, pixel_y)
    cheb_e2 = ChebyshevBoundedField.fit(lsst.geom.Box2I(bbox), x, y, psf_e2, ctrl)
    psf_e2_pts = chebE2.evaluate(pixel_x, pixel_y)

    return psf_size_pts, psf_e1_pts, psf_e2_pts
