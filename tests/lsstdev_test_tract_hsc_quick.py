import unittest
import os
import shutil
import tempfile
import numpy as np

import healsparse
import supreme

import lsst.daf.persistence as dafPersist


class TestTractHscQuickTestCase(unittest.TestCase):
    """
    Tests for running a single tract, quick mode, with HSC data on lsst-dev01
    """
    def test_tract_quick(self):
        """
        Test a single tract in quick mode
        """
        repo = '/datasets/hsc/repo/rerun/RC/w_2020_07/DM-23564'
        tract = 9615
        filter_name = 'HSC-I'

        self.test_dir = tempfile.mkdtemp(dir='./', prefix='TestTractHsc-')

        butler = dafPersist.Butler(repo)
        config = supreme.Configuration(os.path.join('configs/config_quick_tract.yaml'))

        mapper = supreme.TractMapper(butler, config, self.test_dir)
        mapper.run(filter_name, tract)

        # Check that all the patch_inputs are there.  Sum up the expected
        # valid_pixels

        skymap = butler.get('deepCoadd_skyMap')
        tract_info = skymap[tract]
        num_patches = tract_info.getNumPatches()

        em = 'patch_inputs'
        npix = 0
        for ii in range(num_patches[0]):
            for jj in range(num_patches[1]):
                patch = '%d,%d' % (ii, jj)

                map_name = os.path.join(self.test_dir, 'testing_%05d_%s_%s_%s.hs'
                                        % (tract, patch, filter_name, em))
                self.assertTrue(os.path.isfile(map_name))
                m = healsparse.HealSparseMap.read(map_name)

                valid_pixels = m.valid_pixels

                self.assertGreater(valid_pixels.size, 10000)

                npix += valid_pixels.size

        expected_maps = ['exptime_sum']

        for em in expected_maps:
            map_name = os.path.join(self.test_dir, 'testing_%05d_%s_%s.hs'
                                    % (tract, filter_name, em))
            self.assertTrue(os.path.isfile(map_name))
            m = healsparse.HealSparseMap.read(map_name)

            valid_pixels = m.valid_pixels

            self.assertEqual(valid_pixels.size, npix)

            if em == 'exptime_sum':
                self.assertLess(np.max(m.get_values_pix(valid_pixels)), 2201.0)
                self.assertGreater(np.min(m.get_values_pix(valid_pixels)), 199.0)

    def setUp(self):
        self.test_dir = None

    def tearDown(self):
        if self.test_dir is not None:
            if os.path.exists(self.test_dir):
                shutil.rmtree(self.test_dir, True)


if __name__ == '__main__':
    unittest.main()
