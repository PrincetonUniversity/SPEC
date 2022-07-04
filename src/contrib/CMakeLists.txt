list (APPEND spec_sources
	${CMAKE_CURRENT_SOURCE_DIR}/d1mach.f
	${CMAKE_CURRENT_SOURCE_DIR}/dcuhre.f
	${CMAKE_CURRENT_SOURCE_DIR}/i1mach.f
	${CMAKE_CURRENT_SOURCE_DIR}/ilut.f
	${CMAKE_CURRENT_SOURCE_DIR}/iqpack.f
	${CMAKE_CURRENT_SOURCE_DIR}/iters.f
	${CMAKE_CURRENT_SOURCE_DIR}/minpack.f
	${CMAKE_CURRENT_SOURCE_DIR}/rksuite.f
)
set (spec_sources "${spec_sources}" PARENT_SCOPE)
