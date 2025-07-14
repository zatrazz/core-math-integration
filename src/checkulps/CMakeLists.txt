add_executable (checkulps
	        checkulps.cc
		refimpls.cc
		refimpls_mpfr.c
)
set (CMAKE_INCLUDE_CURRENT_DIR on)

find_package(OpenMP REQUIRED)

find_package(PkgConfig REQUIRED)
pkg_check_modules(MPFR REQUIRED mpfr)
find_package(Boost REQUIRED COMPONENTS program_options headers)

add_subdirectory(core-math)

target_include_directories(checkulps PRIVATE "${COMMON_INCLUDE_DIR}")
target_include_directories(checkulps PUBLIC ${MPFR_INCLUDE_DIRS})

target_link_libraries(checkulps PRIVATE Boost::program_options)
target_link_libraries(checkulps PRIVATE ${MPFR_LIBRARIES})
target_link_libraries(checkulps PRIVATE coremathlib)
target_link_libraries(checkulps PRIVATE OpenMP::OpenMP_CXX)
