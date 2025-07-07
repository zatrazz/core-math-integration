CC         = gcc
CXX        = g++
MATH_FLAGS = -fno-finite-math-only -frounding-math -fsignaling-nans

CFLAGS     = -Wall -g -std=gnu11 -D_GNU_SOURCE $(MATH_FLAGS)
CXXFLAGS   = -Wall -g -std=c++23 -D_GNU_SOURCE $(MATH_FLAGS)
ifndef DEBUG
CFLAGS    += -O2 -march=native
CXXFLAGS  += -O2 -march=native
else
CFLAGS    += -O0
CXXFLAGS  += -O0
endif
ifdef ASAN
CFLAGS    += -fsanitize=address
CXXFLAGS  += -fsanitize=address
endif
ifdef UBSAN
CFLAGS    += -fsanitize=undefined
CXXFLAGS  += -fsanitize=undefined
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

checkulps:		checkulps.o \
			refimpls.o \
			refimpls_mpfr.o \
			acos.o \
			acospi.o \
			acosh.o \
			asin.o \
			asinpi.o \
			asinh.o \
			atan.o \
			atan2.o \
			atanh.o \
			atanpi.o \
			cbrt.o \
			cos.o \
			cosh.o \
			cospi.o \
			hypot.o \
			erf.o \
			erfc.o \
			exp.o \
			exp10.o \
			sin.o \
			sinh.o \
			sinpi.o \
			tan.o \
			tanh.o \
			tanpi.o
			$(CXX) $(CXXFLAGS) $(OPENMPFLAGS) -o $@ $^ \
				-lboost_program_options \
				-Wl,-Bstatic -lmpfr -lgmp -Wl,-Bdynamic

atan2.o:		atan2.c
			$(CC) $(CFLAGS) -ffinite-math-only -c -o $@ $^

checkulps.o:		checkulps.cc
			$(CXX) $(CXXFLAGS) $(OPENMPFLAGS) -c -o $@ $^

clean:
			rm -f $(TARGET) *.o
