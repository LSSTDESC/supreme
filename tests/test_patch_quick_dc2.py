import unittest
import os
import tempfile
from collections import OrderedDict

import supreme

import lsst.daf.persistence as dafPersist
import lsst.utils

import supreme_test_base


class PatchQuickDc2TestCase(supreme_test_base.SupremeTestBase):
    """
    Tests for running a single patch, quick mode, with ImSim DC2 data
    """
    @classmethod
    def setUpClass(cls):
        try:
            cls.data_dir = lsst.utils.getPackageDir('supreme_testdata')
        except LookupError:
            raise unittest.SkipTest("supreme_testdata not setup")

        cls.repo = os.path.join(cls.data_dir, 'supreme', 'testdata', 'DC2_test', 'rerun', 'coadd')
        cls.butler = dafPersist.Butler(cls.repo)

    def test_patch_quick(self):
        """
        Test a single RC2 patch in quick mode
        """
        tract = 3828
        filter_name = 'i'
        patch = '2,2'

        self.test_dir = tempfile.mkdtemp(dir='./', prefix='TestPatchImSim-')

        config = supreme.Configuration(os.path.join('configs/config_quick_dc2.yaml'))

        mapper = supreme.MultiMapper(self.butler, config, self.test_dir)
        mapper([tract], [filter_name], [patch], consolidate=False)

        # Check that everything is there
        expected_dict = OrderedDict()
        expected_dict['patch_inputs'] = [2, 1]
        expected_dict['airmass_max'] = [1.0, 1.2, 'float64']
        expected_dict['airmass_min'] = [1.0, 1.2, 'float64']
        expected_dict['airmass_wmean'] = [1.0, 1.2, 'float64']
        expected_dict['background_wmean'] = [2290.0, 2297.0, 'float64']
        expected_dict['bgmean_wmean'] = [2290.0, 2297.0, 'float64']
        expected_dict['exptime_sum'] = [29.0, 61.0, 'float64']
        expected_dict['nexp_sum'] = [0, 3, 'int32']
        expected_dict['psf_size_wmean'] = [1.35, 1.75, 'float64']
        expected_dict['psf_e1_wmean'] = [-0.01, 0.01, 'float64']
        expected_dict['psf_e2_wmean'] = [-0.01, 0.01, 'float64']
        expected_dict['coadd_image_mean'] = [0.0, 2.0, 'float64']
        expected_dict['coadd_variance_mean'] = [0.4, 0.45, 'float64']
        expected_dict['coadd_mask_or'] = [-1, 21000, 'int32']

        self.check_expected_maps_patch(expected_dict, tract, patch, filter_name)


if __name__ == '__main__':
    unittest.main()
