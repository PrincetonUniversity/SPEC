# =============================================================================
# Python wrapper for SPEC using f90wrap and f2py
# =============================================================================

# Load Python and various components
find_package(Python 3.10 REQUIRED
    COMPONENTS Interpreter Development.Module NumPy
)

# Get the path of the `f2py` include directory that will be needed for the library compilation
execute_process(
    COMMAND "${Python_EXECUTABLE}" -c "import numpy.f2py; print(numpy.f2py.get_include())"
    OUTPUT_VARIABLE F2PY_INCLUDE_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
message(STATUS "f2py include dir: ${F2PY_INCLUDE_DIR}")

# Create the list of Fortran files to include in the Python module
# These come from src/CMakeLists.txt via fortran_src_files variable
set(spec4py_f90_src_files ${fortran_src_files})
message(STATUS "Fortran source files for Python wrapper: ${spec4py_f90_src_files}")

# =============================================================================
# Step 1: Preprocess Fortran files using the Fortran compiler preprocessor
# =============================================================================

# Set preprocessor flags for each compiler
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(PP_FLAG "-cpp")
    set(PP_ONLY_FLAG "-E")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(PP_FLAG "-fpp")
    set(PP_ONLY_FLAG "-E")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC|PGI")
    set(PP_FLAG "-Mpreprocess")
    set(PP_ONLY_FLAG "-E")
else()
    message(WARNING "Unknown Fortran compiler: preprocessing may not work")
    set(PP_FLAG "-cpp")
    set(PP_ONLY_FLAG "-E")
endif()

# Preprocess each Fortran source file
foreach(src_file ${spec4py_f90_src_files})
    get_filename_component(basename "${src_file}" NAME_WE)
    set(pp_file "${CMAKE_CURRENT_BINARY_DIR}/${basename}_fpp.F90")
    list(APPEND spec4py_f90_pp_files "${pp_file}")
    add_custom_command(
        COMMAND ${CMAKE_Fortran_COMPILER} ${PP_FLAG} ${PP_ONLY_FLAG} ${src_file}
                -o ${pp_file}
        DEPENDS ${src_file}
        OUTPUT ${pp_file}
        COMMENT "Preprocess ${src_file} -> ${basename}_fpp.F90"
    )
endforeach()

message(STATUS "Preprocessed files: ${spec4py_f90_pp_files}")

# =============================================================================
# Step 2: Run f90wrap to generate Fortran wrapper files
# =============================================================================

# f90wrap generates wrapper files only for modules that have wrappable content.
# These are the known output files based on the SPEC source structure.
# f90wrap_toplevel.f90 contains the standalone subroutines.
set(f90wrap_output_files
    ${CMAKE_CURRENT_BINARY_DIR}/f90wrap_global_m_fpp.f90
    ${CMAKE_CURRENT_BINARY_DIR}/f90wrap_inputlist_m_fpp.f90
    ${CMAKE_CURRENT_BINARY_DIR}/f90wrap_intghs_m_fpp.f90
    ${CMAKE_CURRENT_BINARY_DIR}/f90wrap_msphdf5_m_fpp.f90
    ${CMAKE_CURRENT_BINARY_DIR}/f90wrap_newton_m_fpp.f90
    ${CMAKE_CURRENT_BINARY_DIR}/f90wrap_toplevel.f90
)

# The module name must be "spec_f90wrapped" to be compatible with simsopt
# which imports "spec.spec_f90wrapped"
set(f90wrap_module_name "spec_f90wrapped")
set(f90wrap_python_file "${CMAKE_CURRENT_BINARY_DIR}/${f90wrap_module_name}.py")

add_custom_command(
    COMMAND "${Python_EXECUTABLE}"
            -m f90wrap
            ${spec4py_f90_pp_files}
            --mod-name "${f90wrap_module_name}"
            --kind-map "${CMAKE_CURRENT_SOURCE_DIR}/kind_map"
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    DEPENDS ${spec4py_f90_pp_files}
    OUTPUT ${f90wrap_output_files} ${f90wrap_python_file}
    COMMENT "Generate Fortran interfaces using f90wrap"
    VERBATIM
)

add_custom_target(
    f90wrap_SPEC
    DEPENDS ${f90wrap_output_files} ${f90wrap_python_file}
)

# =============================================================================
# Step 3: Run f2py-f90wrap to generate C interface code
# =============================================================================

set(f2py_module_name "_${f90wrap_module_name}")
set(f2py_c_file "${CMAKE_CURRENT_BINARY_DIR}/${f2py_module_name}module.c")

add_custom_command(
    COMMAND f2py-f90wrap
            -m "${f2py_module_name}"
            --lower
            ${f90wrap_output_files}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    DEPENDS ${f90wrap_output_files}
    OUTPUT ${f2py_c_file}
    COMMENT "Generate C interfaces for Python using f2py-f90wrap"
    VERBATIM
)

add_custom_target(
    f2py_SPEC
    DEPENDS ${f2py_c_file}
)

# =============================================================================
# Step 4: Build the Python extension module
# =============================================================================

# Create the shared library using the generated files
add_library("${f2py_module_name}" MODULE
    ${f2py_c_file}
    "${F2PY_INCLUDE_DIR}/fortranobject.c"
    ${f90wrap_output_files}
)

# Set the library properties for Python extension
set_target_properties("${f2py_module_name}" PROPERTIES
    PREFIX ""
    SUFFIX ".${Python_SOABI}.so"
    C_VISIBILITY_PRESET hidden
    Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules
)

# Tell CMake that the source files are generated
set_source_files_properties(
    ${f2py_c_file}
    ${f90wrap_output_files}
    PROPERTIES GENERATED TRUE
)

target_include_directories("${f2py_module_name}" PRIVATE 
    ${F2PY_INCLUDE_DIR} 
    ${Python_NumPy_INCLUDE_DIRS} 
    ${Python_INCLUDE_DIRS} 
    ${CMAKE_BINARY_DIR}/src
    ${CMAKE_Fortran_MODULE_DIRECTORY}/spec_modules
)

# Compile options for the wrapper Fortran files
target_compile_options("${f2py_module_name}" PRIVATE
    $<$<COMPILE_LANGUAGE:Fortran>:-cpp>
    $<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<Fortran_COMPILER_ID:GNU>>:-ffree-line-length-none>
    $<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<Fortran_COMPILER_ID:GNU>>:-fdefault-real-8>
    $<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<Fortran_COMPILER_ID:GNU>>:-fPIC>
)

target_link_libraries("${f2py_module_name}" PRIVATE 
    Python::NumPy 
    spec
    ${SPEC_LINK_LIB}
)

# Ensure proper build order
add_dependencies("${f2py_module_name}" f90wrap_SPEC f2py_SPEC spec)

# =============================================================================
# Step 5: Installation
# =============================================================================

# The package name is "spec" and the wrapper module inside is "spec_f90wrapped"
# This matches what simsopt expects: "import spec.spec_f90wrapped"
set(spec_package_name "spec")

# Install the existing spec package files from the source directory
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/spec/"
        DESTINATION "${spec_package_name}"
        COMPONENT python_wrapper
        FILES_MATCHING PATTERN "*.py"
)

# Install the generated Python wrapper file (spec_f90wrapped.py)
install(FILES "${f90wrap_python_file}"
        DESTINATION "${spec_package_name}"
        COMPONENT python_wrapper
)

# Install the shared library into the spec package
install(TARGETS "${f2py_module_name}" 
        LIBRARY DESTINATION "${spec_package_name}"
        COMPONENT python_wrapper
)
