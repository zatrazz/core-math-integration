add_executable (checkulps
	        checkulps.cc
		refimpls.cc
		refimpls_binary32_mpfr.c
		refimpls_binary64_mpfr.c
)
set (CMAKE_INCLUDE_CURRENT_DIR on)

find_package(OpenMP REQUIRED)

find_package(PkgConfig REQUIRED)

find_package(Boost REQUIRED COMPONENTS program_options headers)

# We prefere to statically link libmpfr to avoid the overhead of
# TLS access through the default -mtls-variant used by the system
# package.
set(CMAKE_FIND_LIBRARY_SUFFIXES .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
pkg_check_modules(MPFR REQUIRED mpfr)

find_library(GMP_STATIC_LIB libgmp.a REQUIRED)
find_library(MPFR_STATIC_LIB libmpfr.a REQUIRED)

target_include_directories(checkulps PRIVATE "${COMMON_INCLUDE_DIR}")
target_include_directories(checkulps PUBLIC ${MPFR_INCLUDE_DIRS})

target_link_libraries(checkulps PRIVATE Boost::program_options)
target_link_libraries(checkulps PRIVATE OpenMP::OpenMP_CXX)
target_link_libraries(checkulps PRIVATE ${MPFR_STATIC_LIB} ${GMP_STATIC_LIB})
target_link_libraries(checkulps PRIVATE m)

if(APPLE)
    target_link_options(checkulps PRIVATE -undefined dynamic_lookup)
endif()
