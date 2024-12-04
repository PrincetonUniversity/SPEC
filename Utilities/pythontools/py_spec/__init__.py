try:
    from importlib import metadata
except ImportError:
    # Running on pre-3.8 Python; use importlib-metadata package
    import importlib_metadata as metadata

__version__ = metadata.version(__package__ or __name__)

from .ci import test
from .input.spec_namelist import SPECNamelist
from .output.spec import SPECout
from .math.spec_fft import spec_fft
from .math.spec_invfft import spec_invfft
