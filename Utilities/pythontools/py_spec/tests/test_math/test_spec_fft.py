import unittest
import numpy as np
from py_spec.math.spec_fft import spec_fft as fft
from py_spec.math.spec_invfft import spec_invfft as invfft


class fftTests(unittest.TestCase):
    def test_fft_1d(self):
        """
        Test spec_fft on an analytical function, 1D output
        """

        nt = 33 # Number of theta points
        nz = 17 # Number of zeta points
        tarr = np.linspace(0, 2*np.pi, nt, endpoint=False)
        zarr = np.linspace(0, 2*np.pi, nz, endpoint=False)

        tgrid, zgrid = np.meshgrid( tarr, zarr ) # grid

        f = 0.3 \
          + 2*np.cos(tgrid)   \
          - 3*np.cos(2*tgrid-zgrid) \
          + 1*np.cos(2*tgrid+zgrid) \
          - 3.8*np.sin(tgrid+2*zgrid) \
          - 2.4*np.sin(3*zgrid)

        efmn, ofmn, _, _ = fft( tarr, zarr, f, Mpol=2, Ntor=3 )

        places=8
        self.assertAlmostEqual(efmn[ 0], 0.3, places=places)
        self.assertAlmostEqual(efmn[ 7],   2, places=places)
        self.assertAlmostEqual(efmn[13],   1, places=places)
        self.assertAlmostEqual(efmn[15],  -3, places=places)
        self.assertAlmostEqual(ofmn[ 3], 2.4, places=places)
        self.assertAlmostEqual(ofmn[ 5],-3.8, places=places)

    
    def test_fft_2d(self):
        """
        Test spec_fft on an analytical function, 2D output
        """

        nt = 32 # Number of theta points
        nz = 21 # Number of zeta points
        tarr = np.linspace(0, 2*np.pi, nt, endpoint=False)
        zarr = np.linspace(0, 2*np.pi, nz, endpoint=True)

        tgrid, zgrid = np.meshgrid( tarr, zarr ) # grid

        f = 0.3 \
          + 2*np.cos(tgrid)   \
          - 3*np.cos(2*tgrid-zgrid) \
          + 1*np.cos(2*tgrid+zgrid) \
          - 3.8*np.sin(tgrid+2*zgrid) \
          - 2.4*np.sin(3*zgrid)

        Mpol=4
        Ntor=8
        efmn, ofmn = fft( tarr, zarr, f, Mpol=Mpol, Ntor=Ntor, output='2D' )

        places=8
        self.assertAlmostEqual(efmn[Ntor,0], 0.3, places=places)
        self.assertAlmostEqual(efmn[Ntor,1],   2, places=places)
        self.assertAlmostEqual(efmn[Ntor-1,2],   1, places=places)
        self.assertAlmostEqual(efmn[Ntor+1,2],  -3, places=places)
        self.assertAlmostEqual(ofmn[Ntor-2,1],-3.8, places=places)
        self.assertAlmostEqual(ofmn[Ntor+3,0], 2.4, places=places)


    def test_invfft(self):
        """
        Test spec_invfft using analytical function.
        """

        nt = 64 # Number of theta points
        nz = 48 # Number of zeta points
        tarr = np.linspace(0, 2*np.pi, nt, endpoint=False)
        zarr = np.linspace(0, 2*np.pi, nz, endpoint=True)

        tgrid, zgrid = np.meshgrid( tarr, zarr ) # grid

        f = 0.3 \
          + 2*np.cos(tgrid)   \
          - 3*np.cos(2*tgrid-zgrid) \
          + 1*np.cos(2*tgrid+zgrid) \
          - 3.8*np.sin(tgrid+2*zgrid) \
          - 2.4*np.sin(3*zgrid)

        Mpol=8
        Ntor=8
        efmn, ofmn, _im, _in = fft( tarr, zarr, f, Mpol=Mpol, Ntor=Ntor, output='1D' )

        freal = invfft(tarr, zarr, efmn, ofmn, _im, _in)

        self.assertTrue( np.max(np.max(np.abs((freal-f)/f)))<1e-8 )



    