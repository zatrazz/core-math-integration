CC=gcc
CXX=g++
CFLAGS=-O0 -g
CXXFLAGS=-O0 -g

TARGET = \
 acosf_inputs \
 acoshf_inputs \
 log10f_inputs \
 log1pf_inputs \
 erff_inputs \
 erfcf_inputs \
 lgammaf_inputs \
 # TARGET

TESTS = \
 cbrtf_tests \
 # TESTS


all:			$(TARGET) $(TESTS)

acosf_inputs:		acosf_inputs.o acosf.o
			$(CXX) $(CXXFLAGS) -o $@ $^
acosf_inputs.o:		acosf_inputs.cc
			$(CXX) $(CXXFLAGS) -c -o $@ $^
acosf.o:		acosf.c
			$(CC) $(CFLAGS) -c -o $@ $^

acoshf_inputs:		acoshf_inputs.o acoshf.o
			$(CXX) $(CXXFLAGS) -o $@ $^
acoshf_inputs.o:	acoshf_inputs.cc
			$(CXX) $(CXXFLAGS) -c -o $@ $^
acoshf.o:		acoshf.c
			$(CC) $(CFLAGS) -c -o $@ $^


log10f_inputs:		log10f_inputs.o log10f.o powf.o
			$(CXX) $(CXXFLAGS) -o $@ $^
log10f_inputs.o:	log10f_inputs.cc
			$(CXX) $(CXXFLAGS) -c -o $@ $^

log10f.o:		log10f.c
			$(CC) $(CFLAGS) -c -o $@ $^
powf.o:			powf.c
			$(CC) $(CFLAGS) -c -o $@ $^


log1pf_inputs:		log1pf_inputs.o log1pf.o
			$(CXX) $(CXXFLAGS) -o $@ $^
log1pf_inputs.o:	log1pf_inputs.cc
			$(CXX) $(CXXFLAGS) -c -o $@ $^
log1pf.o:		log1pf.c
			$(CC) $(CFLAGS) -c -o $@ $^


cbrtf_tests:		cbrtf_tests.o cbrtf.o
			$(CC) $(CFLAGS) -o $@ $^ -lmpfr -lm
cbrtf_tests.o:		cbrtf_tests.c
			$(CC) $(CFLAGS) -c -o $@ $^
cbrtf.o:		cbrtf.c
			$(CC) $(CFLAGS) -c -o $@ $^


erff_inputs:		erff_inputs.o erff.o
			$(CXX) $(CFLAGS) -o $@ $^
erff_inputs.o:		erff_inputs.cc
			$(CC) $(CFLAGS) -c -o $@ $^
erff.o:			erff.c
			$(CC) $(CFLAGS) -c -o $@ $^


erfcf_inputs:		erfcf_inputs.o erfcf.o
			$(CXX) $(CFLAGS) -o $@ $^
erfcf_inputs.o:		erfcf_inputs.cc
			$(CC) $(CFLAGS) -c -o $@ $^
erfcf.o:			erfcf.c
			$(CC) $(CFLAGS) -c -o $@ $^


lgammaf_inputs:		lgammaf_inputs.o lgammaf.o powf.o
			$(CXX) $(CFLAGS) -o $@ $^
lgammaf_inputs.o:	lgammaf_inputs.cc
			$(CC) $(CFLAGS) -c -o $@ $^
lgammaf.o:		lgammaf.c
			$(CC) $(CFLAGS) -c -o $@ $^
			

clean:
			rm -f $(TARGET) *.o
