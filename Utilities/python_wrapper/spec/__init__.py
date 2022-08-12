import os.path
#print(__file__)
#print(os.path.dirname(__file__))
path_to_spec_f90wrapped = os.path.dirname(__file__)

import sys
if not path_to_spec_f90wrapped in sys.path:
  sys.path.append(path_to_spec_f90wrapped)

#import spec_f90wrapped as spec
