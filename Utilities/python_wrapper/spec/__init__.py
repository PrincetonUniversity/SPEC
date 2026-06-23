import os.path
import sys

# Add current directory to path if not already there
path_to_spec_f90wrapped = os.path.dirname(__file__)
if path_to_spec_f90wrapped not in sys.path:
    sys.path.append(path_to_spec_f90wrapped)

