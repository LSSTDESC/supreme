import os
import numpy as np
import healpy as hp
import healsparse
import esutil

from lsst.pipe.tasks.objectMasks import ObjectMaskCatalog


class RegionMasker(object):
    """
    Make a mask from region files.
    """
    def __init__(self, config):
        """
        Instantiate a RegionMasker

        Parameters
        ----------
        config : `RegionConfiguration`
           supreme region configuration object
        """
        self.config = config

    def __call__(self, maskfile, regionfiles, clobber=False):
        """
        Run the region masker.

        Parameters
        ----------
        maskfile : `str`
           Output mask file
        regionfiles : `list` of `str`
           list of region files to turn into a mask.
        clobber : `bool`, optional
           Clobber any existing maskfile?
        """
        if os.path.isfile(maskfile) and not clobber:
            print("Maskfile %s already exists and clobber is False",
                  maskfile)
            return

        # Read in the region file(s) and convert
        started = False
        for regionfile in regionfiles:
            cat = ObjectMaskCatalog.read(regionfile)

            if not started:
                ctr = 0
                maskcat = np.zeros(len(cat._catalog), dtype=[('ra', 'f8'),
                                                             ('dec', 'f8'),
                                                             ('radius', 'f8'),
                                                             ('width', 'f8'),
                                                             ('height', 'f8')])
                started = True
            else:
                ctr = maskcat.size
                maskcat.resize(maskcat.size + len(cat._catalog), refcheck=False)

            maskcat['ra'][ctr:] = np.rad2deg(cat._catalog['coord_ra'])
            maskcat['dec'][ctr:] = np.rad2deg(cat._catalog['coord_dec'])
            maskcat['radius'][ctr:] = np.rad2deg(cat._catalog['radius'])
            maskcat['width'][ctr:] = np.rad2deg(cat._catalog['width'])
            maskcat['height'][ctr:] = np.rad2deg(cat._catalog['height'])

            # Clear memory
            del cat

        # Split into rendering pixels to conserve memory
        ipnest = hp.ang2pix(self.config.nside_render, maskcat['ra'], maskcat['dec'],
                            lonlat=True, nest=True)
        h, rev = esutil.stat.histogram(ipnest, rev=True)

        gdpix, = np.where(h > 0)
        print('There are %d pixels to render' % (gdpix.size))

        # Figure out approximate coverage to pre-fill
        ipnest_cov = hp.ang2pix(self.config.nside_coverage, maskcat['ra'], maskcat['dec'],
                                lonlat=True, nest=True)
        cov_pixels = np.unique(ipnest_cov)

        # Start the map
        maskmap = healsparse.HealSparseMap.make_empty(self.config.nside_coverage,
                                                      self.config.nside,
                                                      healsparse.WIDE_MASK,
                                                      wide_mask_maxbits=8,
                                                      cov_pixels=cov_pixels)

        # Render all the pixels
        for ctr, i in enumerate(gdpix):
            if (ctr % 100) == 0:
                print('Working on %d' % (ctr))
            i1a = rev[rev[i]: rev[i + 1]]
            self._render_catalog(maskmap, maskcat, i1a)

        # Output
        maskmap.write(maskfile, clobber=clobber)

    def _render_catalog(self, maskmap, maskcat, indices):
        """
        Render part of a maskcat into a mask map.

        Parameters
        ----------
        maskmap : `healsparse.HealSparseMap`
           The map must be of `healsparse.WIDE_MASK` type
        maskcat : `np.ndarray`
           Catalog of region information
        indices : `np.ndarray`
           Array of indices from maskcat to render
        """
        box = np.isnan(maskcat['radius'][indices])
        circ = ~box

        ra = maskcat['ra'][indices[circ]]
        dec = maskcat['dec'][indices[circ]]
        radius = maskcat['radius'][indices[circ]]

        shapes = []
        for i in range(ra.size):
            shapes.append(healsparse.Circle(ra=ra[i],
                                            dec=dec[i],
                                            radius=radius[i],
                                            value=[self.config.default_bit]))

        ra = maskcat['ra'][indices[box]]
        dec = maskcat['dec'][indices[box]]
        width = maskcat['width'][indices[box]] / np.cos(np.deg2rad(maskcat['dec'][indices[box]]))
        height = maskcat['height'][indices[box]]

        for i in range(ra.size):
            shapes.append(healsparse.Polygon(ra=[ra[i] - width[i]/2.,
                                                 ra[i] - width[i]/2.,
                                                 ra[i] + width[i]/2.,
                                                 ra[i] + width[i]/2.],
                                             dec=[dec[i] - height[i]/2.,
                                                  dec[i] + height[i]/2.,
                                                  dec[i] + height[i]/2.,
                                                  dec[i] - height[i]/2.],
                                             value=[self.config.default_bit]))

        healsparse.realize_geom(shapes, maskmap)
