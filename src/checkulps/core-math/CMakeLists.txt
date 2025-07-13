file(GLOB CORE_MATH_SRCS *.c)

# Same as used by CORE-MATH
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fno-finite-math-only -frounding-math -fsignaling-nans")

set_source_files_properties(atan2.c PROPERTIES COMPILE_FLAGS "${CMAKE_C_FLAGS} -ffinite-math-only")

add_library(coremathlib STATIC ${CORE_MATH_SRCS})
