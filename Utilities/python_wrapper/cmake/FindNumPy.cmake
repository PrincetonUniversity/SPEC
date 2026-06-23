# FindNumPy.cmake
# Find NumPy include directories

# Find Python first
find_package(Python COMPONENTS Interpreter NumPy REQUIRED)

# Set NumPy variables for compatibility
set(NumPy_INCLUDE_DIRS ${Python_NumPy_INCLUDE_DIRS})
set(NumPy_FOUND ${Python_NumPy_FOUND})

mark_as_advanced(NumPy_INCLUDE_DIRS)