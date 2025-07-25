include(CheckCCompilerFlag)
include(CheckFunctionExists)

file(GLOB CORE_MATH_SRCS *.c)

check_c_compiler_flag(-fsignaling-nans HAVE_FSIGNALING_NANS)

# Same as used by CORE-MATH
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fno-finite-math-only -frounding-math")
if(HAVE_FSIGNALING_NANS)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsignaling-nans")
endif()

set(CMAKE_REQUIRED_LIBRARIES m)
check_function_exists(roundeven HAVE_ROUNDEVEN)
if (HAVE_ROUNDEVEN)
  list(REMOVE_ITEM CORE_MATH_SRCS "${CMAKE_CURRENT_SOURCE_DIR}/roundeven.c")
endif()
check_function_exists(roundevenf HAVE_ROUNDEVENF)
if (HAVE_ROUNDEVENF)
  list(REMOVE_ITEM CORE_MATH_SRCS "${CMAKE_CURRENT_SOURCE_DIR}/roundevenf.c")
endif()

set_source_files_properties(atan2.c PROPERTIES COMPILE_FLAGS "${CMAKE_C_FLAGS} -ffinite-math-only")

add_library(coremathlib STATIC ${CORE_MATH_SRCS})
