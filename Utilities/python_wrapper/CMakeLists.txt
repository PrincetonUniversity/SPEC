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
message(STATUS "f2py res" ${F2PY_INCLUDE_DIR}
)

# Create the list of Fortran files to include in the Python module
set(spec4py_f90_src_files ${fortran_src_files})

set(f90wrap_output_files ${CMAKE_CURRENT_BINARY_DIR}/f90wrap_global_m.f90
                         ${CMAKE_CURRENT_BINARY_DIR}/f90wrap_inputlist_m.f90
                         ${CMAKE_CURRENT_BINARY_DIR}/f90wrap_intghs_m.f90
                         ${CMAKE_CURRENT_BINARY_DIR}/f90wrap_msphdf5_m.f90
                         ${CMAKE_CURRENT_BINARY_DIR}/f90wrap_newton_m.f90
                         ${CMAKE_CURRENT_BINARY_DIR}/f90wrap_toplevel.f90
)

# Before using `f90wrap`, we need to preprocess the Fortran files because `f90wrap` will otherwise ignore some directives.
# Since there is no way of generalizing the preprocessing, we need to define the correct flags for each compiler.
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
    message(ERROR "Unknown Fortran compiler: preprocessing may not work")
endif()

# TODO: Add the preprocessing options
message(WARNING "Make sure all options are accounted for in the preprocessing step")
foreach(src_file ${spec4py_f90_src_files})
    get_filename_component(basename "${src_file}" NAME_WE)
    set(pp_file "${CMAKE_CURRENT_BINARY_DIR}/${basename}.i90")
    list(APPEND spec4py_f90_pp_files "${pp_file}")
    add_custom_command(
        COMMAND ${CMAKE_Fortran_COMPILER} ${PP_FLAG} ${PP_ONLY_FLAG} ${src_file}
                -o ${pp_file}
        DEPENDS ${src_file}
        OUTPUT ${pp_file}
        COMMENT "Preprocess ${src_file}"
    )
endforeach()

# Create the list of f90wrap output files
# We need to define those filenames to use dependency between the `f90wrap` and `f2py` custom commands.
foreach(file ${spec4py_f90_pp_files})
    get_filename_component(basename "${file}" NAME_WE)
    list(APPEND spec4py_wrap_files "f90wrap_${basename}.f90")
endforeach()

list(TRANSFORM spec4py_wrap_files PREPEND "${CMAKE_CURRENT_BINARY_DIR}/")

add_custom_target(
    f90wrap_SPEC
    DEPENDS ${spec4py_wrap_files}
)

set(f90wrap_module_name "spec4py")
add_custom_command(
    COMMAND "${Python_EXECUTABLE}"
    -m f90wrap
    ${spec4py_f90_pp_files}
    --mod-name "${f90wrap_module_name}"
    --kind-map "${CMAKE_CURRENT_SOURCE_DIR}/kind_map"
    --package
    DEPENDS ${spec4py_f90_pp_files}
    OUTPUT ${spec4py_wrap_files}
    COMMENT "Generate Fortran interfaces using `f90wrap`"
    VERBATIM
)

set(f2py_module_name "_${f90wrap_module_name}")
add_custom_command(
    COMMAND f2py-f90wrap
            -m "${f2py_module_name}"
            --lower
            ${f90wrap_output_files}
    DEPENDS ${spec4py_wrap_files} 
    OUTPUT "${f2py_module_name}module.c"
    COMMENT "Generate C interfaces for Python using `f2py`"
    VERBATIM
)

Python_add_library("${f2py_module_name}" MODULE WITH_SOABI
    "${f2py_module_name}module.c"
    "${F2PY_INCLUDE_DIR}/fortranobject.c"
    "${f90wrap_output_files}"
)

target_include_directories("${f2py_module_name}" PRIVATE ${F2PY_INCLUDE_DIR} ${Python_NumPy_INCLUDE_DIRS} ${Python_INCLUDE_DIRS} ${CMAKE_BINARY_DIR}/src)
target_link_libraries("${f2py_module_name}" PRIVATE Python::NumPy spec_f90)
add_dependencies("${f2py_module_name}" f90wrap_SPEC)

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${f90wrap_module_name}.py" DESTINATION .)
install(TARGETS "${f2py_module_name}" LIBRARY DESTINATION .)
