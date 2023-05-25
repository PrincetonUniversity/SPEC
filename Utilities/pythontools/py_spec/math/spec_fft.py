import numpy as np
import copy
import warnings

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
        If output=='1D', the output is a 4-tuple 1D numpy array, with modes 
            organized as in SPEC (first increasing n, then increasing m). First
            tuple element are the even modes, second are the odd modes, third
            are the poloidal mode number, and fourth the toroidal mode number.
        If output=='2D', the output is tuple of 2D numpy array, with modes (m,n) 
            in element (Ntor+n,m). First tuple elements are the even modes, 
            seconds are the odd modes
    """

    # Check input
    # -----------
    # tarr
    if not isinstance(tarr, np.ndarray):
        raise ValueError('tarr should be a numpy array')
    if tarr.ndim>1:
        raise ValueError('tarr should be one-dimensional')

    nt = tarr.size
    if nt<3:
        raise ValueError('Not enough points in tarr')
    
    # zarr
    if not isinstance(zarr, np.ndarray):
        raise ValueError('zarr should be a numpy array')
    if zarr.ndim>1:
        raise ValueError('zarr should be one-dimensional')
    
    nz = zarr.size
    if nz<3:
        raise ValueError('Not enough points in zarr')
    
    # freal
    if not isinstance(freal, np.ndarray):
        raise ValueError('freal should be a numpy array')
    if freal.ndim!=2:
        raise ValueError('freal should be two-dimensional')
    if freal.shape!=(nz,nt):
        raise ValueError('freal should be of size (nz,nt)')
    
    if tarr[0]+2*np.pi==tarr[-1]:
        warnings.warn('tarr[-1] should not be equal to tarr[0]. Removing last element...')
        tarr = tarr[:-1]
        nt -= 1
        freal = freal[:,:-1]

    if zarr[0]+2*np.pi==zarr[-1]:
        warnings.warn('zarr[-1] should not be equal to zarr[0]. Removing last element...')
        zarr = zarr[:-1]
        nz -= 1
        freal = freal[:-1,:]

    if not np.all(np.abs(np.diff(np.diff(tarr)))<1e-14):
        raise ValueError('Points should be equidistant in tarr')
    if not np.all(np.abs(np.diff(np.diff(zarr)))<1e-14):
        raise ValueError('Points should be equidistant in zarr')

    # Mpol
    M = int(np.floor(nt/2))
    if Mpol is None:
        Mpol = M
    if not isinstance(Mpol, int):
        raise ValueError('Mpol should be an integer')
    if Mpol<1 or Mpol>nt/2:
        raise ValueError('Mpol should be at least 1 and smaller than nt/2')
    
    # Ntor
    N = int(np.floor(nz/2))
    if Ntor is None:
        Ntor = N
    if not isinstance(Ntor, int):
        raise ValueError('Ntor should be an integer')
    if Ntor<1 or Ntor>nz/2:
        raise ValueError('Ntor should be at least 1 and smaller than nt=z/2')
    
    # output
    if not isinstance(output, str):
        raise ValueError('output should be a string')
    if not (output=='1D' or output=='2D'):
        raise ValueError('output should be 1D or 2D')
    
    # Do the fft
    # ----------
    # Actual fft
    fmn = np.fft.fft2( freal ) / (nt*nz)

    # Now construct 2D array. Negative m-modes have to be added to positive 
    # m-modes, i.e. (m,n)+(-m,-n) for even modes, and (m,n)-(-m,-n) for odd
    # modes.
    # Difference if the number of elements (either nt or nz) is odd - this is 
    # due to the internal implementation of numpy.fft.fft.
    efmn = copy.deepcopy(fmn)
    efmn[1:,1:] = fmn[1:,1:] + np.flip(fmn[1:,1:])

    ofmn = copy.deepcopy(fmn)
    ofmn[1:,1:] = fmn[1:,1:] - np.flip(fmn[1:,1:])
    if np.mod(nt,2)==0:
        efmn[0,1:M+1] += np.flip(fmn[0,M:]) #n=0
        ofmn[0,1:M+1] -= np.flip(fmn[0,M:])
    else:
        efmn[0,1:M+1] += np.flip(fmn[0,M+1:]) #n=0
        ofmn[0,1:M+1] -= np.flip(fmn[0,M+1:])

    if np.mod(nz,2)==0:
        efmn[1:N+1,0] += np.flip(fmn[N:,0]) # m=0
        ofmn[1:N+1,0] -= np.flip(fmn[N:,0])
    else:
        efmn[1:N+1,0] += np.flip(fmn[N+1:,0]) # m=0
        ofmn[1:N+1,0] -= np.flip(fmn[N+1:,0])
    
    efmn = (np.absolute(efmn) * np.cos(np.angle(efmn)))[:,0:M+1]
    ofmn = (np.absolute(ofmn) * np.sin(np.angle(ofmn)))[:,0:M+1]

    # Change sign of m>0 odd modes
    ofmn[:,1:] = -ofmn[:,1:]

    # Shift result to be centered around n=0 mode
    efmn = np.fft.fftshift(efmn, axes=0)
    ofmn = np.fft.fftshift(ofmn, axes=0)

    # Invert n -> -n
    if np.mod(nz,2)==0:
        efmn[1:,1:] = np.flip(efmn[1:,1:],axis=0)
        ofmn[1:,1:] = np.flip(ofmn[1:,1:],axis=0)
    else:
        efmn[:,1:] = np.flip(efmn[:,1:],axis=0)
        ofmn[:,1:] = np.flip(ofmn[:,1:],axis=0)


    # Prepare output
    # --------------
    if output=='2D':
        # 2D - truncate efmn, ofmn to the required resolution       
        even_out = efmn[N-Ntor:N+Ntor+1,0:Mpol+1]
        odd_out  = ofmn[N-Ntor:N+Ntor+1,0:Mpol+1]
        return (even_out, odd_out)
    
    elif output=='1D':
        # 1D - construct array mode by mode
        nmn = Ntor+1 + Mpol*(2*Ntor+1)
        even_out = np.zeros((nmn,))
        odd_out = np.zeros((nmn,))
        _im = np.zeros((nmn,), dtype=int)
        _in = np.zeros((nmn,), dtype=int)

        ind = -1
        for mm in range(0,Mpol+1):
            for nn in range(-Ntor,Ntor+1):
                if mm==0 and nn<0:
                    continue
                
                ind += 1

                _im[ind]=mm
                _in[ind]=nn
                even_out[ind] = efmn[N+nn,mm]
                odd_out[ind]  = ofmn[N+nn,mm]

        return (even_out, odd_out, _im, _in)
    
    else:
        raise ValueError('Invalid output')





    
