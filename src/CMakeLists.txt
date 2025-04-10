
# ATM the sources are added as it is. For some files, preprocessing is required. Use custom command to do preprocessing once it is understood.

set(MACROS ${CMAKE_CURRENT_SOURCE_DIR}/macros)

function(preprocess_m4 outvar)
    message(STATUS "preprocess_fortran arguments: ${outvar}, followed by ${ARGN}")
    set(srcs)
    foreach(f ${ARGN})
        message(STATUS "Got file: ${f}")
        # construct output filename
        if(NOT IS_ABSOLUTE "${f}")
            get_filename_component(f "${f}" ABSOLUTE)
        endif()
        file(RELATIVE_PATH r "${CMAKE_CURRENT_SOURCE_DIR}" "${f}")
        get_filename_component(e "${r}" EXT)
        get_filename_component(n "${r}" NAME_WE)
        get_filename_component(p "${r}" PATH)
        string(TOUPPER "${e}" e)
        set(of "${CMAKE_CURRENT_BINARY_DIR}/${n}_m${e}")
        #set(of "${n}_m${e}")
        message(STATUS "Output name: ${of}")
        # preprocess the thing
        add_custom_command(OUTPUT ${of}
            COMMAND m4
            ARGS -P ${MACROS} ${f} >  ${of}
            DEPENDS ${MACROS} ${f}
            COMMENT "Preprocessing ${f}"
            VERBATIM
        )
        list(APPEND srcs "${of}")
    endforeach()
    # return the (preprocessed) sources
    set(${outvar} "${srcs}" PARENT_SCOPE)
endfunction()

set(unprocessed_src_files
    manual.f90
    rzaxis.f90
    packxi.f90
    volume.f90
    coords.f90
    basefn.f90
    memory.f90
    metrix.f90
    ma00aa.f90
    matrix.f90
    spsmat.f90
    spsint.f90
    mp00ac.f90
    ma02aa.f90
    packab.f90
    tr00ab.f90
    curent.f90
    df00ab.f90
    lforce.f90
    intghs.f90
    mtrxhs.f90
    lbpol.f90
    brcast.f90
    dfp100.f90
    dfp200.f90
    dforce.f90
    newton.f90
    casing.f90
    bnorml.f90
    jo00aa.f90
    pp00aa.f90
    pp00ab.f90
    bfield.f90
    stzxyz.f90
    hesian.f90
    ra00aa.f90
    numrec.f90
    preset.f90
    #sphdf5.f90
    global.f90
    inputlist.f90
    #${CMAKE_CURRENT_SOURCE_DIR}/xspech.f90
) # below assumes the .f files are double precision; the CFLAGS = -r8 option is not required;

set(f_src_files
    ${CMAKE_CURRENT_SOURCE_DIR}/dcuhre.f
    ${CMAKE_CURRENT_SOURCE_DIR}/minpack.f
    ${CMAKE_CURRENT_SOURCE_DIR}/iqpack.f
    ${CMAKE_CURRENT_SOURCE_DIR}/rksuite.f
    ${CMAKE_CURRENT_SOURCE_DIR}/i1mach.f
    ${CMAKE_CURRENT_SOURCE_DIR}/d1mach.f
    ${CMAKE_CURRENT_SOURCE_DIR}/ilut.f
    ${CMAKE_CURRENT_SOURCE_DIR}/iters.f
)

preprocess_m4(srcs ${unprocessed_src_files})

set(SPHDF5_FILE "${CMAKE_CURRENT_SOURCE_DIR}/sphdf5.f90")
set(SPHDF5_AWK_FILE "${CMAKE_CURRENT_SOURCE_DIR}/msphdf5.f90")
message(STATUS "sphdf5_FILE is ${SPHDF5_FILE}")
if(AWK MATCHES ".+-NOTFOUND")
    message(FATAL_ERROR "FATAL: awk (and mawk and gawk) could not be found (${AWK}).")
else()
    execute_process(
        COMMAND /bin/sh -c "\"${AWK}\" -v file=sphdf5.f90  \
             '{ gsub(\"__LINE__\", NR); gsub(\"__FILE__\",file); print }'  \
             \"${SPHDF5_FILE}\"" 
        RESULT_VARIABLE AWK_EXITCODE
        OUTPUT_FILE "${SPHDF5_AWK_FILE}"
        )
    message(STATUS "Exit code from awk: ${AWK_EXITCODE}")
endif()
preprocess_m4(sphdf5_src ${SPHDF5_AWK_FILE})
message(STATUS "sphdf5_src is ${sphdf5_src}")

#set(srcs)
#foreach(f ${unprocessed_src_files})
#    # is it a Fortran file?
#    if(f MATCHES "\\.[Ff](9[05])?")
#        message(STATUS "Got fortran file: ${f}")
#        # construct output filename
#        if(NOT IS_ABSOLUTE "${f}")
#            get_filename_component(f "${f}" ABSOLUTE)
#        endif()
#        file(RELATIVE_PATH r "${CMAKE_CURRENT_SOURCE_DIR}" "${f}")
#        get_filename_component(e "${r}" EXT)
#        get_filename_component(n "${r}" NAME_WE)
#        get_filename_component(p "${r}" PATH)
#        #set(of1 "${n}_m${e}")
#        set(of "${CMAKE_CURRENT_BINARY_DIR}/${n}${e}")
#        #set(of "${n}_m${e}")
#        message(STATUS "Output name: ${of}")
#        # preprocess the thing
#        add_custom_command(OUTPUT ${of}
#            COMMAND m4
#            ARGS -P ${MACROS} ${f} >  ${of}
#            DEPENDS ${MACROS} ${f}
#            COMMENT "Preprocessing ${f}"
#            #VERBATIM
#        )
#        list(APPEND srcs "${of}")
#    endif()
#endforeach()


#set_source_files_properties(${m4_output_files} 
#set_source_files_properties(${src} 
#                            PROPERTIES GENERATED true)

#add_custom_target(m4_process ALL
  #DEPENDS ${m4_output_files} 
#  DEPENDS ${src} 
#)
#add_dependencies(spec m4_process)


list(APPEND srcs "${sphdf5_src}")

add_library(spec3p OBJECT ${f_src_files})
target_compile_options(spec3p
    PRIVATE 
        "-cpp"
        $<$<Fortran_COMPILER_ID:GNU>:-ffree-line-length-none>
        # $<$<Fortran_COMPILER_ID:GNU>:-fdefault-real-8>
        $<$<Fortran_COMPILER_ID:GNU>:-fbounds-check>
        # $<$<Fortran_COMPILER_ID:GNU>:-std=legacy>
        $<$<Fortran_COMPILER_ID:GNU>:-fexternal-blas>
        # $<$<Fortran_COMPILER_ID:Intel>:-r8>
)

set(ALLFILES  manual rzaxis packxi volume coords basefn memory metrix ma00aa matrix spsmat spsint mp00ac ma02aa packab tr00ab curent df00ab lforce intghs mtrxhs lbpol brcast dfp100 dfp200 dforce newton casing bnorml jo00aa pp00aa pp00ab bfield stzxyz hesian ra00aa numrec dcuhre minpack iqpack rksuite i1mach d1mach ilut iters sphdf5 preset global xspech)
string(REPLACE ";" " " ALLFILES_STR "${ALLFILES}")

# Build spec executable
set(XSPEC_FILE "${CMAKE_CURRENT_SOURCE_DIR}/xspech.f90")
set(XSPEC_AWK_FILE "${CMAKE_CURRENT_SOURCE_DIR}/mxspech.f90")

if(AWK MATCHES ".+-NOTFOUND")
    message(FATAL_ERROR "FATAL: awk (and mawk and gawk) could not be found (${AWK}).")
else()
    execute_process(
        COMMAND bash -c "date"
        OUTPUT_VARIABLE DATE
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    execute_process(
        COMMAND /bin/sh -c "\"${AWK}\" -v date=\"${DATE}\" -v pwd=\"$ENV{PWD}\" -v macros=\"${MACROS}\" \
            -v fc=\"${CMAKE_Fortran_COMPILER}\" -v flags=\"${COMP_DEFS}\" -v allfiles=\"${ALLFILES_STR}\" \
            'BEGIN{nfiles=split(allfiles,files,\" \")} \
             {if($2==\"COMPILATION\") { \
                   print \"    write(ounit,*)\\\"      :  compiled  : date    = \"date\" ; \\\"\" ; \
                   print \"    write(ounit,*)\\\"      :            : srcdir  = \"pwd\" ; \\\"\" ; \
                   print \"    write(ounit,*)\\\"      :            : macros  = \"macros\" ; \\\"\" ; \
                   print \"    write(ounit,*)\\\"      :            : fc      = \"fc\" ; \\\"\" ; \
                   print \"    write(ounit,*)\\\"      :            : flags   = \"flags\" ; \\\"\" }} \
             {if($2==\"SUMTIME\") {for (i=1;i<=nfiles;i++) print \"   SUMTIME(\"files[i]\")\"}}\
             {if($2==\"PRTTIME\") {for (i=1;i<=nfiles;i++) print \"   PRTTIME(\"files[i]\")\"}}\
             {print}' \"${XSPEC_FILE}\"" 
        RESULT_VARIABLE AWK_EXITCODE
        OUTPUT_FILE "${XSPEC_AWK_FILE}"
        )
    message(STATUS "Exit code from awk: ${AWK_EXITCODE}")
endif()

preprocess_m4(XSPEC_OUT_FILE ${XSPEC_AWK_FILE})
message(STATUS "XSPEC_OUT_FILE is ${XSPEC_OUT_FILE}")

list(APPEND srcs "${XSPEC_OUT_FILE}")
if(SKBUILD)
    set(fortran_src_files "${srcs}" PARENT_SCOPE)
endif()
message(STATUS "srcs variable is ${srcs}")


add_library(spec 
        ${srcs} "$<TARGET_OBJECTS:spec3p>"
)


#target_compile_options(spec PUBLIC "-cpp")
# For gfortran set the -ffree-line-length-none option
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
    target_link_libraries(spec3p PUBLIC OpenMP::OpenMP_Fortran)
    target_compile_definitions(spec3p PUBLIC OPENMP)
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
    ${FFTW3_LIBRARY_TRUE} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES}
)
target_include_directories(spec 
    PUBLIC 
        ${FFTW3_INCLUDE_DIRS} ${LAPACK_INCLUDE_DIRS} ${BLAS_INCLUDE_DIRS}
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

add_executable(xspec ${XSPEC_OUT_FILE})
target_link_libraries(xspec PUBLIC spec)
set_target_properties(xspec PROPERTIES POSITION_INDEPENDENT_CODE ON)

install(TARGETS xspec spec)
