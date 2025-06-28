CC        = gcc
CXX       = g++
CFLAGS    = -g -std=c11 -D_GNU_SOURCE
CXXFLAGS  = -g -std=c++23 -D_GNU_SOURCE
ifndef DEBUG
CFLAGS   += -O2 -march=native
CXXFLAGS += -O2 -march=native
else
CFLAGS   += -O0
CXXFLAGS += -O0
endif

TARGET = \
 randfloatgen \
 checkinputs \
 checkulps \
 # TARGET

ifdef NOOPENMP
OPENMPFLAGS =
else
OPENMPFLAGS = -fopenmp
endif

all:			$(TARGET) $(TESTS)

randfloatgen:		randfloatgen.o
			$(CXX) $(CXXFLAGS) -o $@ $^ -lboost_program_options

checkinputs:		checkinputs.o
			$(CXX) $(CXXFLAGS) -o $@ $^ -lboost_program_options

checkulps:		checkulps.o refimpls.o
			$(CXX) $(CXXFLAGS) $(OPENMPFLAGS) -o $@ $^ -lboost_program_options -lmpfr

checkulps.o:		checkulps.cc
			$(CXX) $(CXXFLAGS) $(OPENMPFLAGS) -c -o $@ $^

clean:
			rm -f $(TARGET) *.o
