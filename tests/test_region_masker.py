import unittest
import os
import glob
import tempfile
import numpy as np
from numpy import testing
import healsparse
import healpy as hp

import supreme
import supreme_test_base


class RegionMaskerTestCase(supreme_test_base.SupremeTestBase):
    """
    Tests for region masker
    """

    def test_single_regfile(self):
        """
        Test a single region file, with box and circle
        """
        self.test_dir = tempfile.mkdtemp(dir='./', prefix='TestRegionMask-')

        config = supreme.RegionConfiguration.load_yaml(os.path.join('./', 'configs',
                                                                    'config_regmask.yaml'))

        masker = supreme.RegionMasker(config)

        regionfiles = sorted(glob.glob(os.path.join('./', 'regions', '*1.reg')))
        outfile = os.path.join(self.test_dir, 'test_reg1.hs')

        masker(outfile, regionfiles)

        mask = healsparse.HealSparseMap.read(outfile)

        pix = mask.valid_pixels

        # Confirm the pixels are what we think they should be
        pix1 = self._get_circ_pix(mask.nside_sparse, 50.0, 25.0, 0.1)
        pix2 = self._get_box_pix(mask.nside_sparse, 50.0, 25.0, 0.05, 0.3)
        test_pixels = np.unique(np.concatenate((pix1, pix2)))

        testing.assert_array_equal(pix, test_pixels)

    def test_two_regfiles(self):
        """
        Test two region files, with box and circle
        """
        self.test_dir = tempfile.mkdtemp(dir='./', prefix='TestRegionMask-')
        config = supreme.RegionConfiguration.load_yaml(os.path.join('./', 'configs',
                                                                    'config_regmask.yaml'))

        masker = supreme.RegionMasker(config)

        regionfiles = sorted(glob.glob(os.path.join('./', 'regions', '*.reg')))
        outfile = os.path.join(self.test_dir, 'test_reg2.hs')

        masker(outfile, regionfiles)

        mask = healsparse.HealSparseMap.read(outfile)

        pix = mask.valid_pixels

        # Confirm the pixels are what we think they should be
        pix1 = self._get_circ_pix(mask.nside_sparse, 50.0, 25.0, 0.1)
        pix2 = self._get_box_pix(mask.nside_sparse, 50.0, 25.0, 0.05, 0.3)
        pix3 = self._get_circ_pix(mask.nside_sparse, 150.0, 25.0, 0.1)
        pix4 = self._get_box_pix(mask.nside_sparse, 150.0, 25.0, 0.05, 0.3)
        test_pixels = np.unique(np.concatenate((pix1, pix2, pix3, pix4)))

        testing.assert_array_equal(pix, test_pixels)

    def _get_circ_pix(self, nside, ra_cent, dec_cent, radius):
        """
        Get healpix pixels overlapping a circle.

        Parameters
        ----------
        nside : `int`
        ra_cent : `float`
        dec_cent : `float`
        radius : `float`

        Returns
        -------
        pixels : `np.ndarray`
        """
        vec = hp.ang2vec(ra_cent, dec_cent, lonlat=True)
        return hp.query_disc(nside, vec, np.deg2rad(radius), nest=True)

    def _get_box_pix(self, nside, ra_cent, dec_cent, width, height):
        """
        Get healpix pixels overlapping a box.

        Parameters
        ----------
        nside : `int`
        ra_cent : `float`
        dec_cent : `float`
        width : `float`
        height : `float`

        Returns
        -------
        pixels : `np.ndarray`
        """
        wid = width / np.cos(np.deg2rad(dec_cent))
        vertices = hp.ang2vec(np.array([ra_cent - wid/2.,
                                        ra_cent - wid/2.,
                                        ra_cent + wid/2.,
                                        ra_cent + wid/2.]),
                              np.array([dec_cent - height/2.,
                                        dec_cent + height/2.,
                                        dec_cent + height/2.,
                                        dec_cent - height/2.]),
                              lonlat=True)
        return hp.query_polygon(nside, vertices, nest=True)


if __name__ == '__main__':
    unittest.main()
