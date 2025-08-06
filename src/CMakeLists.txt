include(CheckIncludeFileCXX)
include(CheckCXXSourceCompiles)

set(COMMON_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include" "${CMAKE_CURRENT_BINARY_DIR}/include")

check_include_file_cxx("print" HAVE_PRINT_HEADER)

check_cxx_source_compiles("
#include <cmath>
int main () {
  int sign;
  volatile double r = lgammaf_r (0.0f, &sign);
}
" HAVE_LGAMMAF_R_PROTOTYPE)

check_cxx_source_compiles("
#include <cmath>
int main () {
  int sign;
  volatile double r = lgamma_r (0.0, &sign);
}
" HAVE_LGAMMA_R_PROTOTYPE)

check_cxx_source_compiles("
#include <expected>
int main() {
    std::expected<int, float> result = 42;
    return result.has_value() ? 0 : 1;
}
" HAVE_STD_EXPECTED)

set(UINT128_TEST_CODE "
int main() {
    __uint128_t value = 0;
    value = 1;
    return 0;
}
")

check_cxx_source_compiles(
  "${UINT128_TEST_CODE}"
  HAVE_UINT128_T)

configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/include/config.h.in"
    "${CMAKE_CURRENT_BINARY_DIR}/include/config.h"
)

add_subdirectory(checkulps)
add_subdirectory(checkinputs)
add_subdirectory(randfloatgen)
