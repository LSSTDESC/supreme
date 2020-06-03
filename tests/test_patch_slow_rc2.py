import unittest
import os
import tempfile
from collections import OrderedDict

import supreme

import lsst.daf.persistence as dafPersist
import lsst.utils

import supreme_test_base


class PatchSlowRc2TestCase(supreme_test_base.SupremeTestBase):
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

    def test_patch_slow(self):
        """
        Test a single RC2 patch in slow mode
        """
        tract = 9697
        filter_name = 'HSC-I'
        patch = '2,2'

        self.test_dir = tempfile.mkdtemp(dir='./', prefix='TestPatchHsc-')

        config = supreme.Configuration.load_yaml(os.path.join('configs/config_slow_rc2.yaml'))

        mapper = supreme.MultiMapper(self.butler, config, self.test_dir)
        mapper([tract], [filter_name], [patch], consolidate=False)

        expected_dict = OrderedDict()
        expected_dict['patch_inputs'] = [2, 1]
        expected_dict['exptime_sum'] = [199.0, 401.0, 'float64']

        self.check_expected_maps_patch(expected_dict, tract, patch, filter_name)


if __name__ == '__main__':
    unittest.main()
