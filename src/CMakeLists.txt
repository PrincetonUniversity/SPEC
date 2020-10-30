
# ATM the sources are added as it is. For some files, preprocessing is required. Use custom command to do preprocessing once it is understood.

set(MACROS ${CMAKE_CURRENT_SOURCE_DIR}/macros)

function(preprocess_m4 outvar)
    message(STATUS "preprocess_fortran arguments: ${outvar}, followed by ${ARGN}")
    set(srcs)
    foreach(f ${ARGN})
        # is it a Fortran file?
        if(f MATCHES "\\.[Ff](9[05])?")
          	message(STATUS "Got fortran file: ${f}")
          	# construct output filename
          	if(NOT IS_ABSOLUTE "${f}")
            	get_filename_component(f "${f}" ABSOLUTE)
       		endif()
        	file(RELATIVE_PATH r "${CMAKE_CURRENT_SOURCE_DIR}" "${f}")
        	get_filename_component(e "${r}" EXT)
        	get_filename_component(n "${r}" NAME_WE)
        	get_filename_component(p "${r}" PATH)
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
     	endif()
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
    sphdf5.f90
    global.f90
    #${CMAKE_CURRENT_SOURCE_DIR}/xspech.f90
) # below assumes the .f files are double precision; the CFLAGS = -r8 option is not required;

set(f_src_fles
    dcuhre.f
    minpack.f
    iqpack.f
    rksuite.f
    i1mach.f
    d1mach.f
    ilut.f
    iters.f
)

#preprocess_m4(m4_output_files ${unprocessed_src_files})

set(srcs)
foreach(f ${unprocessed_src_files})
    # is it a Fortran file?
    if(f MATCHES "\\.[Ff](9[05])?")
        message(STATUS "Got fortran file: ${f}")
        # construct output filename
        if(NOT IS_ABSOLUTE "${f}")
            get_filename_component(f "${f}" ABSOLUTE)
        endif()
        file(RELATIVE_PATH r "${CMAKE_CURRENT_SOURCE_DIR}" "${f}")
        get_filename_component(e "${r}" EXT)
        get_filename_component(n "${r}" NAME_WE)
        get_filename_component(p "${r}" PATH)
        #set(of1 "${n}_m${e}")
        set(of "${CMAKE_CURRENT_BINARY_DIR}/${n}${e}")
        #set(of "${n}_m${e}")
        message(STATUS "Output name: ${of}")
        # preprocess the thing
        add_custom_command(OUTPUT ${of}
            COMMAND m4
            ARGS -P ${MACROS} ${f} >  ${of}
            DEPENDS ${MACROS} ${f}
            COMMENT "Preprocessing ${f}"
            #VERBATIM
        )
        list(APPEND srcs "${of}")
    endif()
endforeach()


#set_source_files_properties(${m4_output_files} 
#set_source_files_properties(${src} 
#                            PROPERTIES GENERATED true)

#add_custom_target(m4_process ALL
  #DEPENDS ${m4_output_files} 
#  DEPENDS ${src} 
#)
#add_dependencies(spec_lib m4_process)

message(STATUS "srcs variable is ${srcs}")

add_library(spec_lib 
        #${m4_output_files}
        ${srcs}
)
target_sources(spec_lib
    PUBLIC
        ${f_src_files}
)
target_compile_options(spec_lib PUBLIC "-cpp")
target_include_directories(spec_lib PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

# Configure OpenMP
target_link_libraries(spec_lib PUBLIC OpenMP)
target_compile_definitions(spec_lib PUBLIC OPENMP)

# Configure HDF5
target_link_libraries(spec_lib PUBLIC ${HDF5_Fortran_LIBRARIES} ${HDF5_C_LIBRARIES})
target_include_directories(spec_lib PUBLIC ${HDF5_C_INCLUDE_DIRS} ${HDF5_Fortran_INCLUDE_DIRS})
