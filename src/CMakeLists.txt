# below assumes the .f files are double precision
# the CFLAGS = -r8 option is not required
set(src_spec_contrib
    ${CMAKE_CURRENT_SOURCE_DIR}/dcuhre.f
    ${CMAKE_CURRENT_SOURCE_DIR}/minpack.f
    ${CMAKE_CURRENT_SOURCE_DIR}/iqpack.f
    ${CMAKE_CURRENT_SOURCE_DIR}/rksuite.f
    ${CMAKE_CURRENT_SOURCE_DIR}/i1mach.f
    ${CMAKE_CURRENT_SOURCE_DIR}/d1mach.f
    ${CMAKE_CURRENT_SOURCE_DIR}/ilut.f
    ${CMAKE_CURRENT_SOURCE_DIR}/iters.f
)

add_library(spec_contrib OBJECT ${src_spec_contrib})
target_compile_options(spec_contrib
    PRIVATE
        #"-cpp"
        #$<$<Fortran_COMPILER_ID:GNU>:-ffree-line-length-none>
        # $<$<Fortran_COMPILER_ID:GNU>:-fdefault-real-8>
        $<$<Fortran_COMPILER_ID:GNU>:-fbounds-check>
        $<$<Fortran_COMPILER_ID:GNU>:-std=legacy>
        $<$<Fortran_COMPILER_ID:GNU>:-fexternal-blas>
        # $<$<Fortran_COMPILER_ID:Intel>:-r8>
)

set(src_spec
    basefn.F90
    bfield.F90
    bnorml.F90
    brcast.F90
    casing.F90
    coords.F90
    curent.F90
    df00ab.F90
    dforce.F90
    dfp100.F90
    dfp200.F90
    global.F90
    h5utils.F90
    hesian.F90
    inputlist.F90
    intghs.F90
    jo00aa.F90
    lbpol.F90
    lforce.F90
    ma00aa.F90
    ma02aa.F90
    manual.F90
    matrix.F90
    memory.F90
    metrix.F90
    mod_kinds.F90
    mp00ac.F90
    mtrxhs.F90
    newton.F90
    numrec.F90
    packab.F90
    packxi.F90
    pp00aa.F90
    pp00ab.F90
    preset.F90
    ra00aa.F90
    rzaxis.F90
    sphdf5.F90
    spsint.F90
    spsmat.F90
    stzxyz.F90
    tr00ab.F90
    volume.F90
    #wa00aa.F90
)

add_library(spec ${src_spec} "$<TARGET_OBJECTS:spec_contrib>")
target_compile_options(spec
    PUBLIC
        "-cpp"
        $<$<Fortran_COMPILER_ID:GNU>:-ffree-line-length-none>
        $<$<Fortran_COMPILER_ID:GNU>:-fdefault-real-8>
        $<$<Fortran_COMPILER_ID:GNU>:-fbounds-check>
        #$<$<Fortran_COMPILER_ID:GNU>:-std=legacy>
        $<$<Fortran_COMPILER_ID:GNU>:-fexternal-blas>
        $<$<Fortran_COMPILER_ID:Intel>:-r8>
)

set_target_properties (spec PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/spec_modules)
target_include_directories(spec PUBLIC ${CMAKE_Fortran_MODULE_DIRECTORY}/spec_modules)

# If using MPI library, link against it
if(MPI_Fortran_FOUND AND NOT SKBUILD)
    target_link_libraries(spec PUBLIC MPI::MPI_Fortran)
endif()

# Add OpenMP
if (OpenMP_Fortran_FOUND AND NOT SKBUILD)
    target_link_libraries(spec PUBLIC OpenMP::OpenMP_Fortran)
    target_compile_definitions(spec PUBLIC OPENMP)
endif()

# Add threads
#target_link_libraries(spec PUBLIC ${CMAKE_THREAD_LIBS_INIT})
# Alternatively
#target_link_libraries(spec PUBLIC Threads::Threads)

# Add libm
# target_link_libraries(spec PUBLIC ${POW_LIBS})

# Add HDF5
#target_link_libraries(spec PUBLIC ${HDF5_Fortran_HL_LIBRARIES} ${HDF5_Fortran_LIBRARIES} ${HDF5_C_LIBRARIES} )
target_link_libraries(spec PUBLIC ${HDF5_Fortran_LIBRARIES} ${HDF5_C_LIBRARIES} )
target_include_directories(spec PUBLIC ${HDF5_C_INCLUDE_DIRS} ${HDF5_Fortran_INCLUDE_DIRS})
#target_include_directories(spec PUBLIC ${HDF5_INCLUDE_DIRS})
#target_link_libraries(spec PUBLIC hdf5::hdf5 hdf5::hdf5_fortran hdf5::hdf5_hl_fortran)

# Add FFTW, LAPACK and BLAS libraries.
# MKL could be used for all the three
target_link_libraries(spec
    PUBLIC
        ${FFTW_LIBRARIES} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES}
)
target_include_directories(spec
    PUBLIC
        ${FFTW_INCLUDE_DIRS} ${LAPACK_INCLUDE_DIRS} ${BLAS_INCLUDE_DIRS}
)

#if(SKBUILD)
#set_target_properties(spec PROPERTIES POSITION_INDEPENDENT_CODE ON)
#endif()

get_property(COMP_DEFS
            TARGET spec
            PROPERTY COMPILE_OPTIONS
)

# Get the spec library linker properties and pass them to python wrapper
get_target_property(SPEC_LINK_LIB spec LINK_LIBRARIES)
message(STATUS "spec linked libraries are ${SPEC_LINK_LIB}")

get_target_property(SPEC_COMPILE_OPTIONS spec COMPILE_OPTIONS)
message(STATUS "spec compile options are ${SPEC_COMPILE_OPTIONS}")

get_target_property(SPEC_COMPILE_DEFS spec COMPILE_DEFINITIONS)
message(STATUS "spec compile definitions are ${SPEC_COMPILE_DEFS}")

# export linker flags for spec to parent scope
# for re-use when building python wrapper
if(SKBUILD)
#    get_target_property(SPEC_LINK_LIB spec LINK_LIBRARIES) # redundant from debug out above
    set(SPEC_LINK_LIB ${SPEC_LINK_LIB} PARENT_SCOPE)
endif()

add_executable(xspec xspech.F90)
target_link_libraries(xspec PUBLIC spec)
set_target_properties(xspec PROPERTIES POSITION_INDEPENDENT_CODE ON)

install(TARGETS xspec spec)
