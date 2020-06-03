import unittest
import os
import tempfile
import yaml

import supreme

import supreme_test_base


class BadConfigTestCase(supreme_test_base.SupremeTestBase):
    """
    Tests for giving bad configs.
    """
    def _write_dict_as_yaml(self, filename, out_dict):
        with open(filename, 'w') as f:
            yaml.dump(out_dict, stream=f)

    @property
    def _minimal_config_dict(self):
        return {'outbase': 'test',
                'map_types': {'exptime': ['sum']}}

    def test_no_mandatory(self):
        """
        Test if we are missing mandatory fields.
        """
        self.test_dir = tempfile.mkdtemp(dir='./', prefix='TestPatchHsc-')
        test_file = os.path.join(self.test_dir, 'test_bad_config.yml')

        config_dict = {}
        self._write_dict_as_yaml(test_file, config_dict)
        self.assertRaises(TypeError, supreme.Configuration.load_yaml, test_file)

        config_dict = {'outbase': 'test'}
        self._write_dict_as_yaml(test_file, config_dict)
        self.assertRaises(TypeError, supreme.Configuration.load_yaml, test_file)

        config_dict = {'map_types': {'exptime': ['sum']}}
        self._write_dict_as_yaml(test_file, config_dict)
        self.assertRaises(TypeError, supreme.Configuration.load_yaml, test_file)

    def test_bad_nside(self):
        """
        Test if we give a bad nside
        """
        self.test_dir = tempfile.mkdtemp(dir='./', prefix='TestPatchHsc-')
        test_file = os.path.join(self.test_dir, 'test_bad_config.yml')

        config_dict = self._minimal_config_dict
        config_dict['nside'] = 0
        self._write_dict_as_yaml(test_file, config_dict)
        self.assertRaises(RuntimeError, supreme.Configuration.load_yaml, test_file)

        config_dict['nside'] = 1025
        self._write_dict_as_yaml(test_file, config_dict)
        self.assertRaises(RuntimeError, supreme.Configuration.load_yaml, test_file)

    def test_bad_map_type(self):
        """
        Test if we give a bad map type
        """
        self.test_dir = tempfile.mkdtemp(dir='./', prefix='TestPatchHsc-')
        test_file = os.path.join(self.test_dir, 'test_bad_config.yml')

        config_dict = self._minimal_config_dict
        config_dict['map_types'] = {'blah': ['sum']}
        self._write_dict_as_yaml(test_file, config_dict)
        self.assertRaises(RuntimeError, supreme.Configuration.load_yaml, test_file)

    def test_bad_map_operation(self):
        """
        Test if we give a bad map operation
        """
        self.test_dir = tempfile.mkdtemp(dir='./', prefix='TestPatchHsc-')
        test_file = os.path.join(self.test_dir, 'test_bad_config.yml')

        config_dict = self._minimal_config_dict
        config_dict['map_types'] = {'exptime': ['blah']}
        self._write_dict_as_yaml(test_file, config_dict)
        self.assertRaises(RuntimeError, supreme.Configuration.load_yaml, test_file)


if __name__ == '__main__':
    unittest.main()
