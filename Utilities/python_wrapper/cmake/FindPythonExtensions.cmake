# FindPythonExtensions.cmake
# Provides functionality for building Python extension modules

# Find Python
find_package(Python COMPONENTS Interpreter Development REQUIRED)

# Get the extension suffix
execute_process(
    COMMAND ${Python_EXECUTABLE} -c "import sysconfig; print(sysconfig.get_config_var('EXT_SUFFIX'))"
    OUTPUT_VARIABLE PYTHON_EXTENSION_MODULE_SUFFIX
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Function to mark a target as a Python extension module
function(python_extension_module target)
    set_target_properties(${target} PROPERTIES
        PREFIX ""
        SUFFIX "${PYTHON_EXTENSION_MODULE_SUFFIX}"
    )
endfunction()

mark_as_advanced(PYTHON_EXTENSION_MODULE_SUFFIX)