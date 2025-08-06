# core-math-integration

Helper tools and programs I use to integrate the [CORE-MATH](https://core-math.gitlabpages.inria.fr/) implementation on glibc:

- **checkinputs**: check the [glibc benchtest input files](https://sourceware.org/git/?p=glibc.git;a=tree;f=benchtests;h=0028936b32775fed8201821028f0a3d350a20506;hb=HEAD) floating-point range.

- **checkulps**: check the accuracy of libm symbol based either on a class of floating-point number (normal or subnormal) or by a random sample in a region.
 
- **randfloatgen**: generate a random floating point number in a specified range in the glibc benchtest input file format.

## Building from source

The project requires a recent C++ compiler that supports C++23. I build and test with gcc/clang from Ubuntu 24 and on macOS using brew llvm (Apple Clang does not support OpenMP).
