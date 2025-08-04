add_executable (checkulps
	        checkulps.cc
		description.cc
		refimpls.cc
		refimpls_binary32_mpfr.c
		refimpls_binary64_mpfr.c
)
set (CMAKE_INCLUDE_CURRENT_DIR on)

find_package(OpenMP REQUIRED)

find_package(PkgConfig REQUIRED)

# For aarch64/x86_64 the default system mpfr usually is built with
# -mtls-dialect=traditional (aarch64) or -mtls-dialect=gnu (x86),
# which incurs in extra overhead for TLS access compared to the
# GNUv2 abi (-mtls-dialect=desc or -mtls-dialect=gnu2).  Linking
# statically avoids the overhead since this only matter for dynamic
# linking
set(CMAKE_FIND_LIBRARY_SUFFIXES .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
pkg_check_modules(MPFR REQUIRED mpfr)
pkg_check_modules(GMP REQUIRED gmp)

find_library(GMP_STATIC_LIB libgmp.a REQUIRED)
find_library(MPFR_STATIC_LIB libmpfr.a REQUIRED)

target_include_directories(checkulps PRIVATE "${COMMON_INCLUDE_DIR}")
target_include_directories(checkulps PUBLIC ${GMP_INCLUDE_DIRS})
target_include_directories(checkulps PUBLIC ${MPFR_INCLUDE_DIRS})

target_link_libraries(checkulps PRIVATE argparse)
target_link_libraries(checkulps PRIVATE nlohmann_json::nlohmann_json)
target_link_libraries(checkulps PRIVATE OpenMP::OpenMP_CXX)
target_link_libraries(checkulps PRIVATE ${MPFR_STATIC_LIB} ${GMP_STATIC_LIB})
target_link_libraries(checkulps PRIVATE m)

if(APPLE)
    target_link_options(checkulps PRIVATE -undefined dynamic_lookup)
endif()
