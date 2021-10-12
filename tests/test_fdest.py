import unittest

import os
import numpy as np
import numpy.testing as testing

import fdest


ROOT = os.path.abspath(os.path.dirname(__file__))


class FdestTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data_dir = 'data'
        cls.atmfile = os.path.join(ROOT,
                                   cls.data_dir,
                                   'des_atmospheres_test.fits')
        cls.ccdfile = os.path.join(ROOT,
                                   cls.data_dir,
                                   'des_y6a1_filter_and_ccd_throughputs_test.fits')

    def test_initialize(self):
        """Test initialization."""
        _ = fdest.FgcmDesTransmission(self.ccdfile, self.atmfile)

    def test_wavelengths(self):
        """Test setting/getting wavelengths."""
        trans = fdest.FgcmDesTransmission(self.ccdfile, self.atmfile)
        wavelengths = trans.get_wavelengths()

        testing.assert_almost_equal(wavelengths[0], 3000.0)
        testing.assert_almost_equal(wavelengths[-1], 11000.0)

        new_wavelengths = np.arange(3000., 11000., 100.0)
        trans.set_wavelengths(new_wavelengths)

        testing.assert_array_almost_equal(trans.get_wavelengths(),
                                          new_wavelengths)

    def test_get_std_transmission(self):
        """Test getting the standard transmission."""
        trans = fdest.FgcmDesTransmission(self.ccdfile, self.atmfile)

        for band in ['g', 'r', 'i', 'z', 'Y']:
            std = trans.get_std_transmission(band)

            testing.assert_almost_equal(std[0], 0.0)
            testing.assert_almost_equal(std[-1], 0.0)
            self.assertGreater(std.max(), 0.5)

    def test_get_transmission(self):
        """Test getting transmission for an exposure/ccd."""
        trans = fdest.FgcmDesTransmission(self.ccdfile, self.atmfile)

        for band in ['g', 'r', 'i', 'z', 'Y']:
            t = trans.get_transmission(band, 226650, 1)

            testing.assert_almost_equal(t[0], 0.0)
            testing.assert_almost_equal(t[-1], 0.0)
            self.assertGreater(t.max(), 0.5)

    def test_get_illegal_transmissions(self):
        """Test getting transmissions for values out of range."""
        trans = fdest.FgcmDesTransmission(self.ccdfile, self.atmfile)

        self.assertRaises(ValueError, trans.get_std_transmission, 'k')
        self.assertRaises(ValueError, trans.get_transmission, 'k', 226650, 1)
        self.assertRaises(ValueError, trans.get_transmission, 'g', 1000000, 1)
        self.assertRaises(ValueError, trans.get_transmission, 'g', 226650, 100)


if __name__ == '__main__':
    unittest.main()
