add_executable (checkinputs
	        checkinputs.cc
)

find_package(Boost REQUIRED COMPONENTS program_options)

target_link_libraries(checkinputs PRIVATE Boost::program_options)
