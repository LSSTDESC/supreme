import unittest
import os
import tempfile
from collections import OrderedDict

import supreme

import lsst.daf.persistence as dafPersist
import lsst.utils

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

    def test_tract_quick(self):
        """
        Test a single RC2 tract in quick mode
        """
        repo = os.path.join(self.data_dir, 'supreme', 'testdata', 'RC2_test', 'rerun', 'coadd')
        tract = 9697
        filter_name = 'HSC-I'

        self.test_dir = tempfile.mkdtemp(dir='./', prefix='TestPatchHsc-')

        butler = dafPersist.Butler(repo)
        config = supreme.Configuration(os.path.join('configs/config_quick_tract_rc2.yaml'))

        mapper = supreme.TractMapper(butler, config, self.test_dir)
        mapper.run(filter_name, tract)

        expected_dict = OrderedDict()
        expected_dict['patch_inputs'] = [4500, '2,1', '2,2']
        expected_dict['exptime_sum'] = [199.0, 401.0, 'float64']

        self.check_expected_maps_tract(expected_dict, tract, filter_name)


if __name__ == '__main__':
    unittest.main()
