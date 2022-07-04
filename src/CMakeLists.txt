add_subdirectory(contrib)
list (APPEND spec_sources
	${CMAKE_CURRENT_SOURCE_DIR}/basefn.F90
	${CMAKE_CURRENT_SOURCE_DIR}/bfield.F90
	${CMAKE_CURRENT_SOURCE_DIR}/bnorml.F90
	${CMAKE_CURRENT_SOURCE_DIR}/brcast.F90
	${CMAKE_CURRENT_SOURCE_DIR}/casing.F90
	${CMAKE_CURRENT_SOURCE_DIR}/coords.F90
	${CMAKE_CURRENT_SOURCE_DIR}/curent.F90
	${CMAKE_CURRENT_SOURCE_DIR}/df00ab.F90
	${CMAKE_CURRENT_SOURCE_DIR}/dforce.F90
	${CMAKE_CURRENT_SOURCE_DIR}/dfp100.F90
	${CMAKE_CURRENT_SOURCE_DIR}/dfp200.F90
	${CMAKE_CURRENT_SOURCE_DIR}/global.F90
	${CMAKE_CURRENT_SOURCE_DIR}/h5utils.F90
	${CMAKE_CURRENT_SOURCE_DIR}/hesian.F90
	${CMAKE_CURRENT_SOURCE_DIR}/inputlist.F90
	${CMAKE_CURRENT_SOURCE_DIR}/intghs.F90
	${CMAKE_CURRENT_SOURCE_DIR}/jo00aa.F90
	${CMAKE_CURRENT_SOURCE_DIR}/lbpol.F90
	${CMAKE_CURRENT_SOURCE_DIR}/lforce.F90
	${CMAKE_CURRENT_SOURCE_DIR}/ma00aa.F90
	${CMAKE_CURRENT_SOURCE_DIR}/ma02aa.F90
	${CMAKE_CURRENT_SOURCE_DIR}/manual.F90
	${CMAKE_CURRENT_SOURCE_DIR}/matrix.F90
	${CMAKE_CURRENT_SOURCE_DIR}/memory.F90
	${CMAKE_CURRENT_SOURCE_DIR}/metrix.F90
	${CMAKE_CURRENT_SOURCE_DIR}/mod_kinds.F90
	${CMAKE_CURRENT_SOURCE_DIR}/mp00ac.F90
	${CMAKE_CURRENT_SOURCE_DIR}/mtrxhs.F90
	${CMAKE_CURRENT_SOURCE_DIR}/newton.F90
	${CMAKE_CURRENT_SOURCE_DIR}/numrec.F90
	${CMAKE_CURRENT_SOURCE_DIR}/packab.F90
	${CMAKE_CURRENT_SOURCE_DIR}/packxi.F90
#	${CMAKE_CURRENT_SOURCE_DIR}/pc00aa.F90 # not de-NAG-ed yet
#	${CMAKE_CURRENT_SOURCE_DIR}/pc00ab.F90 # only used by pc00aa (?)
	${CMAKE_CURRENT_SOURCE_DIR}/pp00aa.F90
	${CMAKE_CURRENT_SOURCE_DIR}/pp00ab.F90
	${CMAKE_CURRENT_SOURCE_DIR}/preset.F90
	${CMAKE_CURRENT_SOURCE_DIR}/ra00aa.F90
	${CMAKE_CURRENT_SOURCE_DIR}/rzaxis.F90
	${CMAKE_CURRENT_SOURCE_DIR}/sphdf5.F90
	${CMAKE_CURRENT_SOURCE_DIR}/spsint.F90
	${CMAKE_CURRENT_SOURCE_DIR}/spsmat.F90
	${CMAKE_CURRENT_SOURCE_DIR}/stzxyz.F90
	${CMAKE_CURRENT_SOURCE_DIR}/tr00ab.F90
	${CMAKE_CURRENT_SOURCE_DIR}/volume.F90
	${CMAKE_CURRENT_SOURCE_DIR}/wa00aa.F90
)
set (spec_sources "${spec_sources}" PARENT_SCOPE)

add_library(spec ${spec_sources})

# For gfortran set the -ffree-line-length-none option
target_compile_options(spec
    PUBLIC
        "-cpp"
        $<$<Fortran_COMPILER_ID:GNU>:-ffree-line-length-none>
#        $<$<Fortran_COMPILER_ID:GNU>:-fdefault-real-8>
#        $<$<Fortran_COMPILER_ID:GNU>:-fbounds-check>
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

# Add HDF5
target_link_libraries(spec PUBLIC ${HDF5_Fortran_LIBRARIES} ${HDF5_C_LIBRARIES} )
target_include_directories(spec PUBLIC ${HDF5_C_INCLUDE_DIRS} ${HDF5_Fortran_INCLUDE_DIRS})

# Add FFTW, LAPACK and BLAS libraries.
# MKL could be used for all the three
target_link_libraries(spec PUBLIC ${FFTW_LIBRARIES} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
target_include_directories(spec PUBLIC ${FFTW_INCLUDE_DIRS} ${LAPACK_INCLUDE_DIRS} ${BLAS_INCLUDE_DIRS})

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
    set(SPEC_LINK_LIB ${SPEC_LINK_LIB} PARENT_SCOPE)
endif()

add_executable(xspec xspech.F90)
target_link_libraries(xspec PUBLIC spec)
set_target_properties(xspec PROPERTIES POSITION_INDEPENDENT_CODE ON)

install(TARGETS xspec spec)
