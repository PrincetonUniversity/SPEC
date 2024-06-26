# import of all SPEC-related python scripts.
__version__ = "3.3.4"

from .ci import test
from .input.spec_namelist import SPECNamelist
from .output.spec import SPECout
from .math.spec_fft import spec_fft
from .math.spec_invfft import spec_invfft
