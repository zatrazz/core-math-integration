//
// Copyright (c) Adhemerval Zanella. All rights reserved.
// Licensed under the MIT license. See LICENSE file in the project root for
// details.
//

#ifndef REFIMPLS_MPFR_H
#define REFIMPLS_MPFR_H

#include "refimpls_modes.h"

#define REFIMPL_F(name)                                                       \
  void ref_##name (double, unsigned, double[REF_NRND]);                       \
  void ref_##name##f (float, unsigned, float[REF_NRND])

#define REFIMPL_FF(name)                                                      \
  void ref_##name (double, double, unsigned, double[REF_NRND]);               \
  void ref_##name##f (float, float, unsigned, float[REF_NRND])

#define REFIMPL_FLLI(name)                                                    \
  void ref_##name (double, long long int, unsigned, double[REF_NRND]);        \
  void ref_##name##f (float, long long int, unsigned, float[REF_NRND])

#define REFIMPL_FPFP(name)                                                    \
  void ref_##name (double, unsigned, double[REF_NRND], double[REF_NRND]);     \
  void ref_##name##f (float, unsigned, float[REF_NRND], float[REF_NRND])

REFIMPL_F (acos);
REFIMPL_F (acosh);
REFIMPL_F (acospi);
REFIMPL_F (asin);
REFIMPL_F (asinh);
REFIMPL_F (asinpi);
REFIMPL_F (atan);
REFIMPL_F (atanh);
REFIMPL_F (atanpi);
REFIMPL_F (cbrt);
REFIMPL_F (cos);
REFIMPL_F (cosh);
REFIMPL_F (cospi);
REFIMPL_F (erf);
REFIMPL_F (erfc);
REFIMPL_F (exp);
REFIMPL_F (exp10);
REFIMPL_F (exp10m1);
REFIMPL_F (exp2);
REFIMPL_F (exp2m1);
REFIMPL_F (expm1);
REFIMPL_F (lgamma);
REFIMPL_F (log);
REFIMPL_F (log1p);
REFIMPL_F (log2);
REFIMPL_F (log2p1);
REFIMPL_F (log10);
REFIMPL_F (log10p1);
REFIMPL_F (rsqrt);
REFIMPL_F (sin);
REFIMPL_F (sinh);
REFIMPL_F (sinpi);
REFIMPL_F (tan);
REFIMPL_F (tanh);
REFIMPL_F (tanpi);
REFIMPL_F (tgamma);

REFIMPL_FF (atan2);
REFIMPL_FF (hypot);
REFIMPL_FF (pow);
REFIMPL_FF (powr);

REFIMPL_FLLI (compoundn);
REFIMPL_FLLI (pown);
REFIMPL_FLLI (rootn);

REFIMPL_FPFP (sincos);

#undef REFIMPL_F
#undef REFIMPL_FF
#undef REFIMPL_FLLI
#undef REFIMPL_FPFP

#endif
