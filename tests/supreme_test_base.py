import unittest
import os
import numpy as np
import numpy.testing as testing
import shutil

import healsparse


class SupremeTestBase(unittest.TestCase):
    def setUp(self):
        self.test_dir = None

    def tearDown(self):
        if self.test_dir is not None:
            if os.path.exists(self.test_dir):
                shutil.rmtree(self.test_dir, True)

    def check_expected_maps_patch(self, expected_dict, tract, patch, filter_name):
        """
        Check for expected maps, ranges, types for a patch.

        Parameters
        ----------
        expected_dict : `OrderedDict`
           Ordered dictionary with keys of map names, values of [min, max].
           The first item must be the patch input with [nccd, width_expected].
        tract : `int`
           Tract number
        patch : `str`
           Patch name
        filter_name : `str`
           Filter name
        """
        mod_times = []
        for em in expected_dict:
            map_name = os.path.join(self.test_dir, '%d' % (tract), 'patches',
                                    'testing_%05d_%s_%s_%s.hs'
                                    % (tract, patch, filter_name, em))
            self.assertTrue(os.path.isfile(map_name))
            mod_times.append(os.path.getmtime(map_name))
            m = healsparse.HealSparseMap.read(map_name)

            valid_pixels = m.valid_pixels
            # The patch_inputs are special, and need to be first
            if em == 'patch_inputs':
                input_valid_pixels = valid_pixels
                metadata = m.metadata

                # Check the metadata
                self.assertTrue(metadata['WIDEMASK'])
                self.assertEqual(metadata['WWIDTH'], expected_dict[em][1])
                nccd = 0
                nvis = 0
                nwt = 0
                nslv = 0
                nssg = 0
                nbgm = 0
                nbgv = 0
                for i in range(10):
                    if 'B%04dCCD' % (i) in metadata:
                        nccd += 1
                    if 'B%04dVIS' % (i) in metadata:
                        nvis += 1
                    if 'B%04dWT' % (i) in metadata:
                        nwt += 1
                    if 'B%04dSLV' % (i) in metadata:
                        nslv += 1
                    if 'B%04dSSG' % (i) in metadata:
                        nssg += 1
                    if 'B%04dBGM' % (i) in metadata:
                        nbgm += 1
                    if 'B%04dBGV' % (i) in metadata:
                        nbgv += 1
                self.assertEqual(nccd, expected_dict[em][0])
                self.assertEqual(nvis, nccd)
                self.assertEqual(nwt, nccd)
                self.assertEqual(nslv, nccd)
                self.assertEqual(nssg, nccd)
                self.assertEqual(nbgm, nccd)
                self.assertEqual(nbgv, nccd)
            else:
                # Make sure we have the same valid pixels as the input map
                testing.assert_array_equal(valid_pixels, input_valid_pixels)

                self.assertEqual(m.dtype.name, expected_dict[em][2])

                self.assertGreater(np.min(m.get_values_pix(valid_pixels)),
                                   expected_dict[em][0])
                self.assertLess(np.max(m.get_values_pix(valid_pixels)),
                                expected_dict[em][1])

        return mod_times

    def check_expected_maps_tract(self, expected_dict, tract, filter_name):
        """
        Check for expected maps, ranges, types for a tract.

        Parameters
        ----------
        expected_dict : `OrderedDict`
           Ordered dictionary with keys of map names, values of [min, max].
           The first item must be the patch input with BLAH.
        tract : `int`
           Tract number
        filter_name : `str`
           Filter name
        """
        mod_times = []
        for em in expected_dict:
            if em == 'patch_inputs':
                for patch in expected_dict[em][1:]:
                    map_name = os.path.join(self.test_dir, '%d' % (tract), 'patches',
                                            'testing_%05d_%s_%s_%s.hs'
                                            % (tract, patch, filter_name, em))

                    self.assertTrue(os.path.isfile(map_name))
                    m = healsparse.HealSparseMap.read(map_name)

                    valid_pixels = m.valid_pixels

                    self.assertGreater(valid_pixels.size, expected_dict[em][0])
            else:
                map_name = os.path.join(self.test_dir, '%d' % (tract),
                                        'testing_%05d_%s_%s.hs'
                                        % (tract, filter_name, em))

                self.assertTrue(os.path.isfile(map_name))
                mod_times.append(os.path.getmtime(map_name))
                m = healsparse.HealSparseMap.read(map_name)

                valid_pixels = m.valid_pixels

                self.assertEqual(m.dtype.name, expected_dict[em][2])

                self.assertGreater(np.min(m.get_values_pix(valid_pixels)),
                                   expected_dict[em][0])
                self.assertLess(np.max(m.get_values_pix(valid_pixels)),
                                expected_dict[em][1])

        return mod_times

    def check_mod_times(self, mod_times1, mod_times2, mode):
        """
        Compare modification times.

        Parameters
        ----------
        mod_times1 : `list`
           First modification times
        mod_times2 : `list`
           Second modification times
        mode : `str`
           'equal' will check times are equal.  'greater` will check
           that mod_times2 > mod_times1
        """
        if mode == 'equal':
            testing.assert_array_almost_equal(mod_times1, mod_times2)
        elif mode == 'greater':
            self.assertTrue(np.all(np.array(mod_times2) > np.array(mod_times1)))
