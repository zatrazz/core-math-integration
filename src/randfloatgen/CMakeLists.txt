add_executable (randfloatgen
	        randfloatgen.cc
)

find_package(Boost REQUIRED COMPONENTS program_options)

target_link_libraries(randfloatgen PRIVATE Boost::program_options)
