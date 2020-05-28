import unittest
import os
import tempfile
from collections import OrderedDict

import supreme

import lsst.daf.persistence as dafPersist
import lsst.utils

import supreme_test_base


class PatchSlowDc2TestCase(supreme_test_base.SupremeTestBase):
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

    def test_patch_slow(self):
        """
        Test a single DC2 patch in slow mode
        """
        tract = 3828
        filter_name = 'i'
        patch = '2,2'

        self.test_dir = tempfile.mkdtemp(dir='./', prefix='TestPatchImSim-')

        config = supreme.Configuration(os.path.join('configs/config_slow_dc2.yaml'))

        mapper = supreme.MultiMapper(self.butler, config, self.test_dir)
        mapper([tract], [filter_name], [patch], consolidate=False)

        expected_dict = OrderedDict()
        expected_dict['patch_inputs'] = [2, 1]
        expected_dict['exptime_sum'] = [29.0, 61.0, 'float64']
        expected_dict['skylevel_wmean'] = [2292.0, 2298.0, 'float64']
        expected_dict['skysigma_wmean'] = [0.2, 0.5, 'float64']

        self.check_expected_maps_patch(expected_dict, tract, patch, filter_name)


if __name__ == '__main__':
    unittest.main()
