add_executable (randfloatgen
	        randfloatgen.cc
)

target_include_directories(randfloatgen PRIVATE "${COMMON_INCLUDE_DIR}")

find_package(Boost REQUIRED COMPONENTS program_options)

target_link_libraries(randfloatgen PRIVATE Boost::program_options)
