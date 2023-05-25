import numpy as np

def spec_fft(tarr: np.array, 
             zarr: np.array,
             freal: np.array,
             Mpol: int=None,
             Ntor: int=None,
            output: str='1D'):
    """
    Fourier transform as in SPEC fortran source

    Args:
     - tarr: 1D numpy array of size nt. Poloidal angle coordinate, each point
             should be equidistant
     - zarr: 1D numpy array of size nz. Toroidal angle coordinate, each point
             should be equidistant
     - freal: 2D numpy array of size (nz,nt). Function, in real space, for which
             Fourier harmonics should be evaluated
     - Mpol: integer. Poloidal resolution. Default is nt/2.
     - Ntor: integer. Toroidal resolution. Default is nz/2.
     - output: either '1D' or '2D'. Determine the form of the output. 

     Output:
        If output=='1D', the output is a 1D numpy array, with modes organized 
            as in SPEC (first increasing n, then increasing m).
        If output=='2D', the output is a 2D numpy array, with modes (m,n) 
            in element (Ntor+n,m)
    """