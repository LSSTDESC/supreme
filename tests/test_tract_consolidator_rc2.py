import unittest
import os
import tempfile
import numpy as np
import healsparse

import supreme
from supreme.utils import op_str_to_code

import supreme_test_base


class TractConsolidateTestCase(supreme_test_base.SupremeTestBase):
    """
    Tests for consolidating tracts, with HSC RC2 config file.
    """
    def test_tract_consolidate_alltracts(self):
        """
        Test consolidating tracts, no explicit specification (all tracts).
        """
        self.test_dir = tempfile.mkdtemp(dir='./', prefix='TestConsolidateHsc-')

        config = supreme.Configuration.load_yaml(os.path.join('./', 'configs',
                                                              'config_consolidate_rc2.yaml'))
        tracts = [100, 200, 500]
        filters = ['HSC-G']
        self._create_fake_maps(config, tracts, filters)

        # Remove one of the airmass_min files (test missing)
        os.remove(os.path.join(self.test_dir, config.tract_relpath(tracts[0]),
                               config.tract_map_filename(filters[0],
                                                         tracts[0],
                                                         'airmass',
                                                         op_str_to_code('min'))))

        # Remove all of the airmass_max files (test all missing)
        for tract in tracts:
            os.remove(os.path.join(self.test_dir, config.tract_relpath(tract),
                                   config.tract_map_filename(filters[0],
                                                             tract,
                                                             'airmass',
                                                             op_str_to_code('max'))))

        # Run the consolidation
        consolidator = supreme.TractConsolidator(config, self.test_dir)
        consolidated_tracts, map_files, map_inputs = consolidator(filters)

        # Make sure the files are there
        nfiles = 0
        for f in filters:
            for i, map_type in enumerate(config.map_types):
                for j, op_str in enumerate(config.map_types[map_type]):
                    op_code = op_str_to_code(op_str)

                    combofile = os.path.join(self.test_dir,
                                             config.consolidated_map_filename(config.outbase,
                                                                              f,
                                                                              map_type,
                                                                              op_code))

                    if map_type == 'airmass' and op_str == 'max':
                        self.assertFalse(os.path.exists(combofile))
                    else:
                        self.assertTrue(os.path.exists(combofile))
                        nfiles += 1

        # Make sure the input/output files are correct
        self.assertEqual(set(tracts), set(consolidated_tracts))
        self.assertEqual(len(map_files), nfiles)
        self.assertEqual(len(map_inputs), nfiles)

        for i in range(len(map_files)):
            if 'airmass_min' in map_files[i]:
                self.assertEqual(len(map_inputs[i]), len(tracts) - 1)
            else:
                self.assertEqual(len(map_inputs[i]), len(tracts))

        # Rerun with clobber=False
        consolidated_tracts, map_files, map_inputs = consolidator(filters, clobber=False)

        # Check that nothing was created.
        self.assertEqual(len(map_files), 0)
        self.assertEqual(len(map_inputs), 0)

        # Rerun with clobber=True
        consolidated_tracts, map_files, map_inputs = consolidator(filters, clobber=True)

        # Check that the input/output files are correct
        self.assertEqual(len(map_files), nfiles)
        self.assertEqual(len(map_inputs), nfiles)

    def test_tract_consolidate_sometracts(self):
        """
        Test consolidating tracts, explicitly specified.
        """
        self.test_dir = tempfile.mkdtemp(dir='./', prefix='TestConsolidateHsc-')

        config = supreme.Configuration.load_yaml(os.path.join('./', 'configs',
                                                              'config_consolidate_rc2.yaml'))
        tracts = [100, 200, 500]
        run_tracts = [100, 200]
        filters = ['HSC-G']
        self._create_fake_maps(config, tracts, filters)

        # Run the consolidation, test override outputbase
        consolidator = supreme.TractConsolidator(config, self.test_dir)
        consolidated_tracts, map_files, map_inputs = consolidator(filters, tracts=run_tracts,
                                                                  outputbase='sometracts')

        # Make sure the files are there
        nfiles = 0
        for f in filters:
            for i, map_type in enumerate(config.map_types):
                for j, op_str in enumerate(config.map_types[map_type]):
                    op_code = op_str_to_code(op_str)

                    combofile = os.path.join(self.test_dir,
                                             config.consolidated_map_filename('sometracts',
                                                                              f,
                                                                              map_type,
                                                                              op_code))

                    self.assertTrue(os.path.exists(combofile))
                    nfiles += 1

        # Make sure the input/output files are correct
        self.assertEqual(set(run_tracts), set(consolidated_tracts))
        self.assertEqual(len(map_files), nfiles)
        self.assertEqual(len(map_inputs), nfiles)

        for i in range(len(map_files)):
            self.assertEqual(len(map_inputs[i]), len(run_tracts))

    def _create_fake_maps(self, config, tracts, filters):
        """
        Create fake maps

        Parameters
        ----------
        config : `supreme.Configuration`
        tracts : `list` of `int`
        filters : `list` of `str`
        """
        for tract in tracts:
            tract_path = os.path.join(self.test_dir, config.tract_relpath(tract))
            os.makedirs(tract_path)

            for f in filters:
                for i, map_type in enumerate(config.map_types):
                    for j, op_str in enumerate(config.map_types[map_type]):
                        op_code = op_str_to_code(op_str)

                        fname = os.path.join(tract_path,
                                             config.tract_map_filename(f,
                                                                       tract,
                                                                       map_type,
                                                                       op_code))
                        if map_type == 'nexp':
                            dtype = np.int32
                            value = 1
                        else:
                            dtype = np.float64
                            value = 1.0
                        fake_map = healsparse.HealSparseMap.make_empty(nside_coverage=32,
                                                                       nside_sparse=4096,
                                                                       dtype=dtype)
                        fake_map[i] = value
                        fake_map.write(fname)


if __name__ == '__main__':
    unittest.main()
