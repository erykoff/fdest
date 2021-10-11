import numpy as np
import fitsio
import scipy.interpolate as interpolate


class FgcmDesTransmission(object):
    """Class to return S_obs(lambda) for DES FGCM tables.

    Parameters
    ----------
    ccd_file : `str`
        File with ccd throughput data.
    atm_file : `str`
        File with atmosphere throughput data.
    """
    def __init__(self, ccd_file, atm_file):
        self._ccd_data = self._read_ccd_file(ccd_file)
        self._atm_data = self._read_atm_file(atm_file)

        self._wavelengths = None
        self._std_ifunc = None

    def set_wavelengths(self, wavelengths):
        """Set the wavelengths to return for the transmission.

        This will cache interpolation.

        Parameters
        ----------
        wavelengths : `np.ndarray`
            Wavelength array (Angstroms)
        """
        if self._wavelengths is not None:
            if len(wavelengths) == len(self._wavelengths):
                if np.allclose(wavelengths == self._wavelengths):
                    # This is a match, we don't need to do anything
                    return

        self._wavelengths = np.atleast_1d(wavelengths)

        # Do ccd interpolation here.
        self._band_tput_interp = {}
        self._band_ccd_interp = {}
        for band in self.bands:
            ifunc = interpolate.interp1d(self._ccd_data[band]['lambda'],
                                         self._ccd_data[band]['throughput_avg'])
            self._band_tput_interp[band] = np.clip(ifunc(self._wavelengths), 0.0, 1e100)
            self._band_ccd_interp[band] = np.zeros((self._wavelengths.size, self.nccd))
            for j in range(self.nccd):
                ifunc = interpolate.interp1d(self._ccd_data[band]['lambda'],
                                             self._ccd_data[band]['throughput_ccd'][:, j])
                self._band_ccd_interp[band][:, j] = np.clip(ifunc(self._wavelengths, j))

    def get_wavelengths(self):
        """Get the wavelengths in the cache.

        Returns
        -------
        wavelengths : `np.ndarray`
            Wavelength array (Angstroms)
        """
        return self._wavelengths

    def get_transmission(self, band, expnum, ccdnum, wavelengths=None):
        """Get the transmission S(lambda).

        Return the transmission for a given band, exposure, and ccd.

        Parameters
        ----------
        band : `str`
            Band to get transmission.
        expnum : `int`
            Exposure number to get transmission.
        ccdnum : `int`
            CCD number to get transmission.
        wavelengths : `np.ndarray`, optional
            Set new wavelengths and interpolation; otherwise
            use cached values.  Units are Angstroms.

        Returns
        -------
        transmission : `np.ndarray`
            Transmission as a function of wavelength.
        """
        if wavelengths is not None:
            self.set_wavelengths(wavelengths)

        u, = np.where(self._atm_data['expnum'] == expnum)
        if len(u) == 0:
            raise ValueError(f"Exposure {expnum} not found in atm table.")

        ccd_index = ccdnum - 1
        if ccd_index < 0 or ccd_index >= self.nccd:
            raise ValueError(f"ccdnum {ccdnum} out of range.")

        if band not in self._band_ccd_interp:
            raise ValueError(f"band {band} not in throughput table.")

        ifunc = interpolate.interp1d(self._atm_wavelengths,
                                     self._atm_data['throughput'][u, :])
        atm = np.clip(ifunc(self._wavelengths), 0.0, 1e100)

        return atm*self._band_ccd_interp[band][:, ccd_index]

    def get_std_transmission(self, band, wavelengths=None):
        """Get the standard transmission S(lambda).

        Return the standard transmission for a given band.

        Parameters
        ----------
        band : `str`
            Band to get standard transmission
        wavelengths : `np.ndarray`, optional
            Set new wavelengths and interpolation; otherwise
            use cached values.  Units are Angstroms.

        Returns
        -------
        std_transmission : `np.ndarray`
            Standard transmission as a function of wavelength.
        """
        if self._std_ifunc is None:
            self._std_ifunc = interpolate.interp1d(self._atm_wavelengths,
                                                   self._atm_std)
        if wavelengths is not None:
            self.set_wavelengths(wavelengths)

        if band not in self._band_tput_interp:
            raise ValueError(f"band {band} not in throughput table.")

        atm = np.clip(self._std_ifunc(self._wavelengths), 0.0, 1e100)

        return atm*self._band_tput_interp[band]

    def _read_ccd_file(self, ccd_file):
        """Read the ccd file.

        Parameters
        ----------
        ccd_file : `str`
            File with ccd throughput data.
        """
        self._ccd_data = {}
        self.bands = []
        self.nccd = None
        with fitsio.FITS(ccd_file) as fits:
            for hdu in fits:
                parts = hdu.get_extname().split('_')

                self._ccd_data[parts[0]] = hdu.read()
                self.bands.append(parts[0])
                if self.nccd is None:
                    self.nccd = self._ccd_data[parts[0]]['throughput_ccd'].shape[1]

    def _read_atm_file(self, atm_file):
        """Read the atmosphere file.

        Parameters
        ----------
        atm_file : `str`
            File with atmosphere throughput data.
        """
        atm_data = fitsio.read(atm_file)
        # The 0th row has the wavelengths
        self._atm_wavelengths = atm_data['throughput'][0, :]
        # The 1st row has the standard atmosphere
        self._atm_std = atm_data['throughput'][1, :]
        # And the rest of the rows are per-exposure
        self._atm_data = atm_data[2:, :]
