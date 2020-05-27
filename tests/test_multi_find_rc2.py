import unittest
import os

import supreme

import lsst.daf.persistence as dafPersist
import lsst.utils

import supreme_test_base


class MultiRc2TestCase(supreme_test_base.SupremeTestBase):
    """
    Tests for running a multi-mode, with HSC RC2 data
    """
    @classmethod
    def setUpClass(cls):
        try:
            cls.data_dir = lsst.utils.getPackageDir('supreme_testdata')
        except LookupError:
            raise unittest.SkipTest("supreme_testdata not setup")

        cls.repo = os.path.join(cls.data_dir, 'supreme', 'testdata', 'RC2_test', 'rerun', 'coadd')
        cls.butler = dafPersist.Butler(cls.repo)

    def test_find_patch(self):
        """
        Find a single patch.
        """
        tract = 9697
        filter_name = 'HSC-I'

        mapper = supreme.MultiMapper(self.butler, None, './')
        multi_dict = mapper.run([tract], [filter_name], patches=['2,2', '5,5'], find_only=True)

        # Should only return the one patch, '2,2'
        self.assertEqual(len(multi_dict), 1)
        self.assertTrue(tract in multi_dict)
        self.assertEqual(len(multi_dict[tract]), 1)
        self.assertTrue('HSC-I' in multi_dict[tract])
        self.assertEqual(len(multi_dict[tract]['HSC-I']), 1)
        self.assertEqual(multi_dict[tract]['HSC-I'][0], '2,2')

    def test_find_filter(self):
        """
        Find a filter.
        """
        tract = 9697
        filter_names = ['HSC-G', 'HSC-R', 'HSC-I', 'HSC-Z', 'HSC-Y']

        mapper = supreme.MultiMapper(self.butler, None, './')
        multi_dict = mapper.run([tract], filter_names, find_only=True)

        # Should only return the one filter, 'HSC-I', with 2 patches
        self.assertEqual(len(multi_dict), 1)
        self.assertTrue(tract in multi_dict)
        self.assertEqual(len(multi_dict[tract]), 1)
        self.assertTrue('HSC-I' in multi_dict[tract])
        self.assertEqual(len(multi_dict[tract]['HSC-I']), 2)


if __name__ == '__main__':
    unittest.main()
