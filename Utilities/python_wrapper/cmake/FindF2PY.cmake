# FindF2PY.cmake
# Find f2py executable and include directories

# Find Python first
find_package(Python COMPONENTS Interpreter Development NumPy REQUIRED)

# Find f2py executable
execute_process(
    COMMAND ${Python_EXECUTABLE} -c "import numpy.f2py; print(numpy.f2py.get_include())"
    OUTPUT_VARIABLE F2PY_INCLUDE_DIRS
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
)

# Find f2py executable path
execute_process(
    COMMAND ${Python_EXECUTABLE} -c "import sys; import numpy.f2py as f2py; print(f2py.main.__file__.replace('__main__.py', 'f2py.py') if hasattr(f2py, 'main') else sys.executable + ' -m numpy.f2py')"
    OUTPUT_VARIABLE F2PY_EXECUTABLE_PATH
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
)

# Set F2PY_EXECUTABLE - use python -m numpy.f2py for reliability
set(F2PY_EXECUTABLE "${Python_EXECUTABLE}" "-m" "numpy.f2py")

# Set include directories
if(NOT F2PY_INCLUDE_DIRS)
    set(F2PY_INCLUDE_DIRS ${Python_NumPy_INCLUDE_DIRS})
endif()

# Mark as found if we have Python and NumPy
if(Python_FOUND AND Python_NumPy_FOUND)
    set(F2PY_FOUND TRUE)
else()
    set(F2PY_FOUND FALSE)
endif()

mark_as_advanced(F2PY_EXECUTABLE F2PY_INCLUDE_DIRS)