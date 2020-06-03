import unittest
import os
import tempfile
from collections import OrderedDict

import supreme

import lsst.daf.persistence as dafPersist
import lsst.utils

import supreme_test_base


class PatchQuickRc2TestCase(supreme_test_base.SupremeTestBase):
    """
    Tests for running a single patch, quick mode, with HSC RC2 data
    """
    @classmethod
    def setUpClass(cls):
        try:
            cls.data_dir = lsst.utils.getPackageDir('supreme_testdata')
        except LookupError:
            raise unittest.SkipTest("supreme_testdata not setup")

        cls.repo = os.path.join(cls.data_dir, 'supreme', 'testdata', 'RC2_test', 'rerun', 'coadd')
        cls.butler = dafPersist.Butler(cls.repo)

    def test_patch_quick(self):
        """
        Test a single RC2 patch in quick mode
        """
        tract = 9697
        filter_name = 'HSC-I'
        patch = '2,2'

        self.test_dir = tempfile.mkdtemp(dir='./', prefix='TestPatchHsc-')

        config = supreme.Configuration.load_yaml(os.path.join('configs/config_quick_rc2.yaml'))

        mapper = supreme.MultiMapper(self.butler, config, self.test_dir)
        mapper([tract], [filter_name], [patch], consolidate=False)

        # Check that everything is there
        expected_dict = OrderedDict()
        expected_dict['patch_inputs'] = [2, 1]
        expected_dict['airmass_max'] = [1.15, 1.3, 'float64']
        expected_dict['airmass_min'] = [1.15, 1.3, 'float64']
        expected_dict['airmass_wmean'] = [1.15, 1.3, 'float64']
        expected_dict['background_wmean'] = [690.0, 790.0, 'float64']
        expected_dict['bgmean_wmean'] = [690.0, 790.0, 'float64']
        expected_dict['exptime_sum'] = [199.0, 401.0, 'float64']
        expected_dict['nexp_sum'] = [0, 3, 'int32']
        expected_dict['psf_size_wmean'] = [1.15, 1.45, 'float64']
        expected_dict['psf_e1_wmean'] = [0.0, 0.15, 'float64']
        expected_dict['psf_e2_wmean'] = [-0.15, 0.15, 'float64']
        expected_dict['skylevel_wmean'] = [690.0, 790.0, 'float64']
        expected_dict['skysigma_wmean'] = [8.4, 9.0, 'float64']
        expected_dict['coadd_image_mean'] = [-0.1, 15.0, 'float64']
        expected_dict['coadd_variance_mean'] = [0.005, 0.03, 'float64']
        expected_dict['coadd_mask_or'] = [-1, 54000, 'int32']

        self.check_expected_maps_patch(expected_dict, tract, patch, filter_name)


if __name__ == '__main__':
    unittest.main()
