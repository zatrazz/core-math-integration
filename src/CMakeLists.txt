include(CheckIncludeFileCXX)
include(CheckCXXSourceCompiles)

set(COMMON_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include" "${CMAKE_CURRENT_BINARY_DIR}/include")

check_include_file_cxx("print" HAVE_PRINT_HEADER)

check_cxx_source_compiles("
#include <cmath>
int main () {
  volatile double r = lgammaf_r (0.0);
}
" HAVE_LGAMMAF_R_PROTOTYE)

check_cxx_source_compiles("
#include <cmath>
int main () {
  volatile double r = lgamma_r (0.0f);
}

" HAVE_LGAMMA_R_PROTOTYE)
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/include/config.h.in"
    "${CMAKE_CURRENT_BINARY_DIR}/include/config.h"
)

add_subdirectory(checkulps)
add_subdirectory(checkinputs)
add_subdirectory(randfloatgen)
