import unittest
import os
import shutil
import tempfile
import numpy.testing as testing
import numpy as np

import healsparse
import supreme

import lsst.daf.persistence as dafPersist


class TestPatchDC2QuickTestCase(unittest.TestCase):
    """
    Tests for running a single patch, quick mode, with DC2 data on nersc
    """
    def test_patch_quick(self):
        """
        Test a single patch in quick mode
        """
        repo = '/global/cfs/cdirs/lsst/production/DC2_ImSim/Run2.2i/desc_dm_drp/' \
            'v19.0.0-v1/rerun/run2.2i-coadd-wfd-dr3-v1'
        tract = 2907
        filter_name = 'i'
        patch = '0,0'

        self.test_dir = tempfile.mkdtemp(dir='./', prefix='TestPatchDC2-')

        butler = dafPersist.Butler(repo)
        config = supreme.Configuration(os.path.join('configs/config_quick_dc2.yaml'))

        mapper = supreme.PatchMapper(butler, config, self.test_dir)
        mapper.run(filter_name, tract, patch, return_values_list=False)

        # Check that everything is there
        expected_maps = ['patch_inputs',
                         'airmass_max', 'airmass_min', 'airmass_wmean',
                         'background_wmean', 'bgmean_wmean',
                         'exptime_sum', 'nexp_sum',
                         'psf_size_wmean', 'psf_e1_wmean', 'psf_e2_wmean',
                         'coadd_image_mean', 'coadd_variance_mean', 'coadd_mask_or']

        for em in expected_maps:
            map_name = os.path.join(self.test_dir, 'testing_%05d_%s_%s_%s.hs'
                                    % (tract, patch, filter_name, em))
            self.assertTrue(os.path.isfile(map_name))
            m = healsparse.HealSparseMap.read(map_name)

            valid_pixels = m.valid_pixels
            if em == 'patch_inputs':
                input_valid_pixels = valid_pixels
                metadata = m.metadata

                # Check the metadata
                self.assertTrue(metadata['WIDEMASK'])
                self.assertEqual(metadata['WWIDTH'], 8)
                nccd = 0
                nvis = 0
                nwt = 0
                nslv = 0
                nssg = 0
                nbgm = 0
                nbgv = 0
                for i in range(100):
                    if 'B%04dCCD' % (i) in metadata:
                        nccd += 1
                    if 'B%04dVIS' % (i) in metadata:
                        nvis += 1
                    if 'B%04dWT' % (i) in metadata:
                        nwt += 1
                    if 'B%04dSLV' % (i) in metadata:
                        nslv += 1
                    if 'B%04dSSG' % (i) in metadata:
                        nssg += 1
                    if 'B%04dBGM' % (i) in metadata:
                        nbgm += 1
                    if 'B%04dBGV' % (i) in metadata:
                        nbgv += 1

                self.assertEqual(nccd, 62)
                self.assertEqual(nvis, nccd)
                self.assertEqual(nwt, nccd)
                self.assertEqual(nslv, nccd)
                self.assertEqual(nssg, nccd)
                self.assertEqual(nbgm, nccd)
                self.assertEqual(nbgv, nccd)
            else:
                testing.assert_array_equal(valid_pixels, input_valid_pixels)

            if em == 'airmass_max' or em == 'airmass_min' or em == 'armass_wmean':
                self.assertLess(np.max(m.get_values_pix(valid_pixels)), 1.5)
                self.assertGreater(np.max(m.get_values_pix(valid_pixels)), 1.0)
            elif em == 'background_wmean' or em == 'bgmean_wmean' or em == 'skylevel_wmean':
                self.assertLess(np.max(m.get_values_pix(valid_pixels)), 3200.0)
                self.assertGreater(np.min(m.get_values_pix(valid_pixels)), 1000.0)
            elif em == 'exptime_sum':
                self.assertLess(np.max(m.get_values_pix(valid_pixels)), 841.0)
                self.assertGreater(np.min(m.get_values_pix(valid_pixels)), 29.0)
            elif em == 'nexp_sum':
                self.assertEqual(m.dtype.name, 'int32')
                self.assertLess(np.max(m.get_values_pix(valid_pixels)), 28.1)
                self.assertGreater(np.min(m.get_values_pix(valid_pixels)), 0.9)
            elif em == 'psf_size_wmean':
                self.assertLess(np.max(m.get_values_pix(valid_pixels)), 2.3)
                self.assertGreater(np.min(m.get_values_pix(valid_pixels)), 1.5)
            elif em == 'psf_e1_wmean' or em == 'psf_e2_wmean':
                self.assertLess(np.max(m.get_values_pix(valid_pixels)), 0.01)
                self.assertGreater(np.min(m.get_values_pix(valid_pixels)), -0.025)
            elif em == 'coadd_image_mean':
                self.assertLess(np.max(m.get_values_pix(valid_pixels)), 60.0)
                self.assertGreater(np.min(m.get_values_pix(valid_pixels)), -1.0)
            elif em == 'coadd_variance_mean':
                self.assertLess(np.max(m.get_values_pix(valid_pixels)), 1.0)
                self.assertGreater(np.min(m.get_values_pix(valid_pixels)), 0.01)
            elif em == 'coadd_mask_or':
                self.assertEqual(m.dtype.name, 'int32')
                self.assertLess(np.max(m.get_values_pix(valid_pixels)), 56000)
                self.assertGreater(np.min(m.get_values_pix(valid_pixels)), -1)

    def setUp(self):
        self.test_dir = None

    def tearDown(self):
        if self.test_dir is not None:
            if os.path.exists(self.test_dir):
                shutil.rmtree(self.test_dir, True)


if __name__ == '__main__':
    unittest.main()
