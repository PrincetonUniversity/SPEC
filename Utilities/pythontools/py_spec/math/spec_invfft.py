import numpy as np

def spec_invfft(tarr:np.ndarray, zarr:np.ndarray, 
                efmn:np.ndarray, ofmn:np.ndarray, 
                _im:np.ndarray, _in:np.ndarray):
    """
    Inverse Fourier transform as in SPEC Fortran source.

    Args:
     - tarr: 1D numpy array of size nt. Poloidal angle coordinate.
     - zarr: 1D numpy array of size nz. Toroidal angle coordinate.
     - efmn: 1D numpy array of size nmn. Even mode numbers.
     - ofmn: 1D numpy array of size nmn. Odd mode numbers.
     - _im:  1D numpy array of size nmn. Poloidal mode numbers
     - _in:  1D numpy array of size nmn. Toroidal mode numbers (multiples of Nfp)

    Output:
     - freal: 2D numpy array of size nz x nt. Function evaluation in real space.
    """

    # Check input
    # -----------
    # tarr
    if not isinstance(tarr, np.ndarray):
        raise ValueError('tarr should be a numpy array')
    if tarr.ndim>1:
        raise ValueError('tarr should be one-dimensional')

    nt = tarr.size

    # zarr
    if not isinstance(zarr, np.ndarray):
        raise ValueError('zarr should be a numpy array')
    if zarr.ndim>1:
        raise ValueError('zarr should be one-dimensional')
    
    nz = zarr.size

    # efmn
    if not isinstance(efmn, np.ndarray):
        raise ValueError('efmn should be a numpy array')
    if efmn.ndim>1:
        raise ValueError('efmn should be one-dimensional')
    
    nmn = efmn.size

    # ofmn
    if not isinstance(ofmn, np.ndarray):
        raise ValueError('ofmn should be a numpy array')
    if ofmn.ndim>1:
        raise ValueError('ofmn should be one-dimensional')
    if ofmn.size!=nmn:
        raise ValueError('ofmn should have the same size as efmn')
    
    # _im
    if not isinstance(_im, np.ndarray):
        raise ValueError('_im should be a numpy array')
    if _im.ndim>1:
        raise ValueError('_im should be one-dimensional')
    if _im.size!=nmn:
        raise ValueError('_im should have the same size as efmn')
    if not np.all(_im>=0):
        raise ValueError('_im should not have any negative values')
    
    # _in
    if not isinstance(_in, np.ndarray):
        raise ValueError('_in should be a numpy array')
    if _in.ndim>1:
        raise ValueError('_in should be one-dimensional')
    if _in.size!=nmn:
        raise ValueError('_in should have the same size as efmn')
        
    # Inverse Fourier transform
    # -------------------------

    # Generate grid:
    tgrid, zgrid = np.meshgrid(tarr, zarr)

    # Evaluate function:
    freal = np.zeros((nz,nt))
    for mm, nn, emode, omode in zip(_im, _in, efmn, ofmn):
        freal += emode * np.cos(mm*tgrid - nn*zgrid) \
               + omode * np.sin(mm*tgrid - nn*zgrid)
    
    # Output
    return freal



