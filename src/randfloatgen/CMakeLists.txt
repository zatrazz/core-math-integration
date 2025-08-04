add_executable (randfloatgen
	        randfloatgen.cc
)

target_include_directories(randfloatgen PRIVATE "${COMMON_INCLUDE_DIR}")

target_link_libraries(randfloatgen PRIVATE argparse)
