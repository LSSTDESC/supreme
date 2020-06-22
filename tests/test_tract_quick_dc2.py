import unittest
import os
import tempfile
from collections import OrderedDict

import supreme

import lsst.daf.persistence as dafPersist
import lsst.utils
import lsst.obs.base

import supreme_test_base


class TractQuickRc2TestCase(supreme_test_base.SupremeTestBase):
    """
    Tests for running a single tract, quick mode, with HSC RC2 data
    """
    @classmethod
    def setUpClass(cls):
        try:
            cls.data_dir = lsst.utils.getPackageDir('supreme_testdata')
        except LookupError:
            raise unittest.SkipTest("supreme_testdata not setup")

        lsst.obs.base.FilterDefinitionCollection.reset()

        cls.repo = os.path.join(cls.data_dir, 'supreme', 'testdata', 'DC2_test', 'rerun', 'coadd')
        cls.butler = dafPersist.Butler(cls.repo)

    def test_tract_quick(self):
        """
        Test a single RC2 tract in quick mode
        """
        tract = 3828
        filter_name = 'i'

        self.test_dir = tempfile.mkdtemp(dir='./', prefix='TestPatchImSim-')

        config = supreme.Configuration.load_yaml(os.path.join('configs/config_quick_tract_dc2.yaml'))

        mapper = supreme.MultiMapper(self.butler, config, self.test_dir, ncores=2)
        mapper([tract], [filter_name])

        expected_dict = OrderedDict()
        expected_dict['patch_inputs'] = [14000, '2,2']
        expected_dict['exptime_sum'] = [29.0, 61.0, 'float64']

        self.check_expected_maps_tract(expected_dict, tract, filter_name)

    def test_tract_quick_nofilter(self):
        """
        Test building a tract that doesn't exist
        """
        tract = 3828
        filter_name = 'r'

        self.test_dir = tempfile.mkdtemp(dir='./', prefix='TestPatchImSim-')

        config = supreme.Configuration.load_yaml(os.path.join('configs/config_quick_tract_dc2.yaml'))

        mapper = supreme.MultiMapper(self.butler, config, self.test_dir, ncores=1)
        mapper([tract], [filter_name])


if __name__ == '__main__':
    unittest.main()
