CC=gcc
CXX=g++
ifndef DEBUG
CFLAGS=-O2 -g -std=c11 -D_GNU_SOURCE
CXXFLAGS=-O2 -g -std=c++20
else
CFLAGS=-O0 -g -std=c11 -D_GNU_SOURCE
CXXFLAGS=-O0 -g -std=c++20
endif

TARGET = \
 randfloatgen \
 checkinputs \
 # TARGET

all:			$(TARGET) $(TESTS)

randfloatgen:		randfloatgen.o
			$(CXX) $(CXXFLAGS) -o $@ $^ -lboost_program_options

checkinputs:		checkinputs.o
			$(CXX) $(CXXFLAGS) -o $@ $^ -lboost_program_options

clean:
			rm -f $(TARGET) *.o
