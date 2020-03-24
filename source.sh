#
# source.sh
# This file sets the PATH to the current directory, which 
# is the directory where the xspec file is located. 
# 
# It also updates the $PYTHONPATH variable so that the 
# regression testing can be executed. 
#

# Test if this is the folder where xspec lives then export paths or exit with error
if ( [ -e ./xspec ] && \
     [ -e ./xspech.f90 ] ); then 
  export PATH=${PATH}:$PWD
  export PYTHONPATH=${PYTHONPATH}:$PWD/Utilities/pythontools
else 
  echo "xspec executable does not exist, will not modify your path"
  exit 1
fi


