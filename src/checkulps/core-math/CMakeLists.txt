include(CheckCCompilerFlag)

file(GLOB CORE_MATH_SRCS *.c)

check_c_compiler_flag(-fsignaling-nans HAVE_FSIGNALING_NANS)

# Same as used by CORE-MATH
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fno-finite-math-only -frounding-math")
if(HAVE_FSIGNALING_NANS)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsignaling-nans")
endif()

set_source_files_properties(atan2.c PROPERTIES COMPILE_FLAGS "${CMAKE_C_FLAGS} -ffinite-math-only")

add_library(coremathlib STATIC ${CORE_MATH_SRCS})
