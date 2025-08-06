add_executable (checkinputs
	        checkinputs.cc
)

target_include_directories(checkinputs PRIVATE "${COMMON_INCLUDE_DIR}")

target_link_libraries(checkinputs PRIVATE argparse)
