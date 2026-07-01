//
// Copyright (c) Adhemerval Zanella. All rights reserved.
// Licensed under the MIT license. See LICENSE file in the project root for
// details.
//

#include <float.h>
#include <math.h>
#include <stdint.h>
// NB: stdint should be included prior mpfr.h to define intmax_t, which
// mpfr.h uses to define mpfr_pown
#include <mpfr.h>

#include "refimpls_modes.h"

// The reference implementation only expect finite numbers, so no NaN/Inf
// handling is required.
//
// Each function computes the transcendental once at INTERNAL_PRECISION (two
// bits more than the target format) in round-toward-zero, then rounds that
// single high-precision result to double for every rounding mode selected in
// MASK.  MPFR's round-toward-zero result plus its exact ternary is folded into
// a sticky bit (round-to-odd) by set_sticky(); rounding a round-to-odd value
// to the narrower double format in any mode is double-rounding safe, so the
// per-mode mpfr_get_d() calls yield the same values as computing each mode
// independently (and mpfr_get_d() applies the double exponent range, including
// gradual underflow to subnormals and overflow to infinity).

enum
{
  INTERNAL_PRECISION = DBL_MANT_DIG + 2
};

static const mpfr_rnd_t ref_rnd_modes[REF_NRND]
    = { MPFR_RNDN, MPFR_RNDU, MPFR_RNDD, MPFR_RNDZ };

// Per-thread scratch reused across calls to avoid a malloc/free of the MPFR
// limbs (and the sticky-bit mpz) on every evaluation.  Initialized lazily on
// first use and intentionally never freed (the process owns it until exit).
static __thread struct
{
  int inited;
  mpfr_t a, b, c;
  mpz_t sticky;
} scr;

static void
scratch_init (void)
{
  if (__builtin_expect (!scr.inited, 0))
    {
      mpfr_init2 (scr.a, INTERNAL_PRECISION);
      mpfr_init2 (scr.b, INTERNAL_PRECISION);
      mpfr_init2 (scr.c, INTERNAL_PRECISION);
      mpz_init (scr.sticky);
      scr.inited = 1;
    }
}

// Fold the round-toward-zero ternary INEX into a sticky bit (round-to-odd):
// if the result was inexact, force the least significant bit of the mantissa.
static void
set_sticky (mpfr_t r, int inex)
{
  if (inex == 0 || !mpfr_number_p (r))
    return;
  mpfr_exp_t e = mpfr_get_z_2exp (scr.sticky, r);
  if (mpz_sgn (scr.sticky) < 0)
    {
      mpz_neg (scr.sticky, scr.sticky);
      mpz_setbit (scr.sticky, 0);
      mpz_neg (scr.sticky, scr.sticky);
    }
  else
    mpz_setbit (scr.sticky, 0);
  mpfr_set_z_2exp (r, scr.sticky, e, MPFR_RNDN);
}

// Round the high-precision value HI (computed round-toward-zero, INEX its
// ternary) to double for every rounding mode selected in MASK.
static void
round_all (mpfr_t hi, int inex, unsigned mask, double out[REF_NRND])
{
  set_sticky (hi, inex);
  for (int i = 0; i < REF_NRND; i++)
    if (mask & (1u << i))
      out[i] = mpfr_get_d (hi, ref_rnd_modes[i]);
}

// Store the mode-independent value V into every slot selected in MASK.  Used
// for special cases (NaN/Inf/zero) whose result does not depend on rounding.
static void
broadcast (double v, unsigned mask, double out[REF_NRND])
{
  for (int i = 0; i < REF_NRND; i++)
    if (mask & (1u << i))
      out[i] = v;
}

void
ref_acos (double x, unsigned mask, double out[REF_NRND])
{
  scratch_init ();
  mpfr_set_d (scr.a, x, MPFR_RNDN);
  int inex = mpfr_acos (scr.a, scr.a, MPFR_RNDZ);
  round_all (scr.a, inex, mask, out);
}

void
ref_acosh (double x, unsigned mask, double out[REF_NRND])
{
  union
  {
    uint64_t u;
    double f;
  } ix = { .f = x };
  if (__builtin_expect (ix.u <= 0x3ff0000000000000ull, 0))
    {
      if (ix.u == 0x3ff0000000000000ull)
	{
	  broadcast (0, mask, out);
	  return;
	}
      broadcast (__builtin_nan ("x<1"), mask, out);
      return;
    }
  if (__builtin_expect (ix.u >= 0x7ff0000000000000ull, 0))
    {
      uint64_t aix = ix.u << 1;
      if (ix.u == 0x7ff0000000000000ull
	  || aix > ((uint64_t) 0x7ff << DBL_MANT_DIG))
	{
	  broadcast (x, mask, out); // +inf or nan
	  return;
	}
      broadcast (__builtin_nan ("x<1"), mask, out);
      return;
    }

  scratch_init ();
  mpfr_set_d (scr.a, x, MPFR_RNDN);
  int inex = mpfr_acosh (scr.a, scr.a, MPFR_RNDZ);
  round_all (scr.a, inex, mask, out);
}

void
ref_acospi (double x, unsigned mask, double out[REF_NRND])
{
  scratch_init ();
  mpfr_set_d (scr.a, x, MPFR_RNDN);
  int inex = mpfr_acospi (scr.a, scr.a, MPFR_RNDZ);
  round_all (scr.a, inex, mask, out);
}

void
ref_asin (double x, unsigned mask, double out[REF_NRND])
{
  scratch_init ();
  mpfr_set_d (scr.a, x, MPFR_RNDN);
  int inex = mpfr_asin (scr.a, scr.a, MPFR_RNDZ);
  round_all (scr.a, inex, mask, out);
}

void
ref_asinh (double x, unsigned mask, double out[REF_NRND])
{
  scratch_init ();
  mpfr_set_d (scr.a, x, MPFR_RNDN);
  int inex = mpfr_asinh (scr.a, scr.a, MPFR_RNDZ);
  round_all (scr.a, inex, mask, out);
}

void
ref_asinpi (double x, unsigned mask, double out[REF_NRND])
{
  scratch_init ();
  mpfr_set_d (scr.a, x, MPFR_RNDN);
  int inex = mpfr_asinpi (scr.a, scr.a, MPFR_RNDZ);
  round_all (scr.a, inex, mask, out);
}

void
ref_atan (double x, unsigned mask, double out[REF_NRND])
{
  scratch_init ();
  mpfr_set_d (scr.a, x, MPFR_RNDN);
  int inex = mpfr_atan (scr.a, scr.a, MPFR_RNDZ);
  round_all (scr.a, inex, mask, out);
}

void
ref_atan2 (double y, double x, unsigned mask, double out[REF_NRND])
{
  scratch_init ();
  mpfr_set_d (scr.b, x, MPFR_RNDN);
  mpfr_set_d (scr.c, y, MPFR_RNDN);
  int inex = mpfr_atan2 (scr.a, scr.c, scr.b, MPFR_RNDZ);
  round_all (scr.a, inex, mask, out);
}

void
ref_atanh (double x, unsigned mask, double out[REF_NRND])
{
  scratch_init ();
  mpfr_set_d (scr.a, x, MPFR_RNDN);
  int inex = mpfr_atanh (scr.a, scr.a, MPFR_RNDZ);
  round_all (scr.a, inex, mask, out);
}

void
ref_atanpi (double x, unsigned mask, double out[REF_NRND])
{
  scratch_init ();
  mpfr_set_d (scr.a, x, MPFR_RNDN);
  int inex = mpfr_atanpi (scr.a, scr.a, MPFR_RNDZ);
  round_all (scr.a, inex, mask, out);
}

void
ref_cbrt (double x, unsigned mask, double out[REF_NRND])
{
  scratch_init ();
  mpfr_set_d (scr.a, x, MPFR_RNDN);
  int inex = mpfr_cbrt (scr.a, scr.a, MPFR_RNDZ);
  round_all (scr.a, inex, mask, out);
}

void
ref_compoundn (double x, long long int y, unsigned mask, double out[REF_NRND])
{
  scratch_init ();
  mpfr_set_d (scr.b, x, MPFR_RNDN);
  int inex = mpfr_compound_si (scr.a, scr.b, y, MPFR_RNDZ);
  round_all (scr.a, inex, mask, out);
}

void
ref_cos (double x, unsigned mask, double out[REF_NRND])
{
  scratch_init ();
  mpfr_set_d (scr.a, x, MPFR_RNDN);
  int inex = mpfr_cos (scr.a, scr.a, MPFR_RNDZ);
  round_all (scr.a, inex, mask, out);
}

void
ref_cosh (double x, unsigned mask, double out[REF_NRND])
{
  scratch_init ();
  mpfr_set_d (scr.a, x, MPFR_RNDN);
  int inex = mpfr_cosh (scr.a, scr.a, MPFR_RNDZ);
  round_all (scr.a, inex, mask, out);
}

void
ref_cospi (double x, unsigned mask, double out[REF_NRND])
{
  if (isnan (x))
    {
      broadcast (x, mask, out);
      return;
    }
  if (isinf (x))
    {
      broadcast (__builtin_nan (""), mask, out);
      return;
    }
  scratch_init ();
  mpfr_set_d (scr.a, x, MPFR_RNDN);
  int inex = mpfr_cospi (scr.a, scr.a, MPFR_RNDZ);
  round_all (scr.a, inex, mask, out);
}

void
ref_erf (double x, unsigned mask, double out[REF_NRND])
{
  scratch_init ();
  mpfr_set_d (scr.a, x, MPFR_RNDN);
  int inex = mpfr_erf (scr.a, scr.a, MPFR_RNDZ);
  round_all (scr.a, inex, mask, out);
}

void
ref_erfc (double x, unsigned mask, double out[REF_NRND])
{
  scratch_init ();
  mpfr_set_d (scr.a, x, MPFR_RNDN);
  int inex = mpfr_erfc (scr.a, scr.a, MPFR_RNDZ);
  round_all (scr.a, inex, mask, out);
}

void
ref_exp10 (double x, unsigned mask, double out[REF_NRND])
{
  scratch_init ();
  mpfr_set_d (scr.a, x, MPFR_RNDN);
  int inex = mpfr_exp10 (scr.a, scr.a, MPFR_RNDZ);
  round_all (scr.a, inex, mask, out);
}

void
ref_exp10m1 (double x, unsigned mask, double out[REF_NRND])
{
  scratch_init ();
  mpfr_set_d (scr.a, x, MPFR_RNDN);
  int inex = mpfr_exp10m1 (scr.a, scr.a, MPFR_RNDZ);
  round_all (scr.a, inex, mask, out);
}

void
ref_exp2 (double x, unsigned mask, double out[REF_NRND])
{
  scratch_init ();
  mpfr_set_d (scr.a, x, MPFR_RNDN);
  int inex = mpfr_exp2 (scr.a, scr.a, MPFR_RNDZ);
  round_all (scr.a, inex, mask, out);
}

void
ref_exp2m1 (double x, unsigned mask, double out[REF_NRND])
{
  scratch_init ();
  mpfr_set_d (scr.a, x, MPFR_RNDN);
  int inex = mpfr_exp2m1 (scr.a, scr.a, MPFR_RNDZ);
  round_all (scr.a, inex, mask, out);
}

void
ref_exp (double x, unsigned mask, double out[REF_NRND])
{
  scratch_init ();
  mpfr_set_d (scr.a, x, MPFR_RNDN);
  int inex = mpfr_exp (scr.a, scr.a, MPFR_RNDZ);
  round_all (scr.a, inex, mask, out);
}

void
ref_expm1 (double x, unsigned mask, double out[REF_NRND])
{
  scratch_init ();
  mpfr_set_d (scr.a, x, MPFR_RNDN);
  int inex = mpfr_expm1 (scr.a, scr.a, MPFR_RNDZ);
  round_all (scr.a, inex, mask, out);
}

void
ref_hypot (double x, double y, unsigned mask, double out[REF_NRND])
{
  /* since MPFR does not distinguish between quiet/signaling NaN,
     we have to deal with them separately to apply the IEEE rules */
  union
  {
    double f;
    uint64_t u;
  } xi = { .f = x }, yi = { .f = y };
  if ((xi.u << 1) < (0xfffull << 52)
      && (xi.u << 1) > (0x7ffull << DBL_MANT_DIG)) // x = sNAN
    {
      broadcast (x + y, mask, out); // will return qNAN
      return;
    }
  if ((yi.u << 1) < (0xfffull << 52)
      && (yi.u << 1) > (0x7ffull << DBL_MANT_DIG)) // y = sNAN
    {
      broadcast (x + y, mask, out); // will return qNAN
      return;
    }
  if ((xi.u << 1) == 0)
    { // x = +/-0
      yi.u = (yi.u << 1) >> 1;
      broadcast (yi.f, mask, out);
      return;
    }
  if ((yi.u << 1) == 0)
    { // y = +/-0
      xi.u = (xi.u << 1) >> 1;
      broadcast (xi.f, mask, out);
      return;
    }

  scratch_init ();
  mpfr_set_d (scr.b, x, MPFR_RNDN);
  mpfr_set_d (scr.c, y, MPFR_RNDN);
  int inex = mpfr_hypot (scr.a, scr.b, scr.c, MPFR_RNDZ);
  round_all (scr.a, inex, mask, out);
}

void
ref_lgamma (double x, unsigned mask, double out[REF_NRND])
{
  int sign;
  scratch_init ();
  mpfr_set_d (scr.a, x, MPFR_RNDN);
  int inex = mpfr_lgamma (scr.a, &sign, scr.a, MPFR_RNDZ);
  round_all (scr.a, inex, mask, out);
}

void
ref_log (double x, unsigned mask, double out[REF_NRND])
{
  scratch_init ();
  mpfr_set_d (scr.a, x, MPFR_RNDN);
  int inex = mpfr_log (scr.a, scr.a, MPFR_RNDZ);
  round_all (scr.a, inex, mask, out);
}

void
ref_log1p (double x, unsigned mask, double out[REF_NRND])
{
  scratch_init ();
  mpfr_set_d (scr.a, x, MPFR_RNDN);
  int inex = mpfr_log1p (scr.a, scr.a, MPFR_RNDZ);
  round_all (scr.a, inex, mask, out);
}

void
ref_log2 (double x, unsigned mask, double out[REF_NRND])
{
  scratch_init ();
  mpfr_set_d (scr.a, x, MPFR_RNDN);
  int inex = mpfr_log2 (scr.a, scr.a, MPFR_RNDZ);
  round_all (scr.a, inex, mask, out);
}

void
ref_log2p1 (double x, unsigned mask, double out[REF_NRND])
{
  scratch_init ();
  mpfr_set_d (scr.a, x, MPFR_RNDN);
  int inex = mpfr_log2p1 (scr.a, scr.a, MPFR_RNDZ);
  round_all (scr.a, inex, mask, out);
}

void
ref_log10 (double x, unsigned mask, double out[REF_NRND])
{
  scratch_init ();
  mpfr_set_d (scr.a, x, MPFR_RNDN);
  int inex = mpfr_log10 (scr.a, scr.a, MPFR_RNDZ);
  round_all (scr.a, inex, mask, out);
}

void
ref_log10p1 (double x, unsigned mask, double out[REF_NRND])
{
  scratch_init ();
  mpfr_set_d (scr.a, x, MPFR_RNDN);
  int inex = mpfr_log10p1 (scr.a, scr.a, MPFR_RNDZ);
  round_all (scr.a, inex, mask, out);
}

void
ref_rsqrt (double x, unsigned mask, double out[REF_NRND])
{
  /* mpfr_rec_sqrt differs from IEEE 754-2019: IEEE 754-2019 says that
     rsqrt(-0) should give -Inf, whereas mpfr_rec_sqrt(-0) gives +Inf */
  if (x == 0.0 && 1.0 / x < 0)
    {
      broadcast (1.0 / x, mask, out);
      return;
    }
  scratch_init ();
  mpfr_set_d (scr.a, x, MPFR_RNDN);
  int inex = mpfr_rec_sqrt (scr.a, scr.a, MPFR_RNDZ);
  round_all (scr.a, inex, mask, out);
}

void
ref_sin (double x, unsigned mask, double out[REF_NRND])
{
  scratch_init ();
  mpfr_set_d (scr.a, x, MPFR_RNDN);
  int inex = mpfr_sin (scr.a, scr.a, MPFR_RNDZ);
  round_all (scr.a, inex, mask, out);
}

void
ref_sincos (double x, unsigned mask, double sinp[REF_NRND],
	    double cosp[REF_NRND])
{
  scratch_init ();
  mpfr_set_d (scr.a, x, MPFR_RNDN);
  int inex_sin = mpfr_sin (scr.a, scr.a, MPFR_RNDZ);
  round_all (scr.a, inex_sin, mask, sinp);

  mpfr_set_d (scr.b, x, MPFR_RNDN);
  int inex_cos = mpfr_cos (scr.b, scr.b, MPFR_RNDZ);
  round_all (scr.b, inex_cos, mask, cosp);
}

void
ref_sinh (double x, unsigned mask, double out[REF_NRND])
{
  scratch_init ();
  mpfr_set_d (scr.a, x, MPFR_RNDN);
  int inex = mpfr_sinh (scr.a, scr.a, MPFR_RNDZ);
  round_all (scr.a, inex, mask, out);
}

void
ref_sinpi (double x, unsigned mask, double out[REF_NRND])
{
  scratch_init ();
  mpfr_set_d (scr.a, x, MPFR_RNDN);
  int inex = mpfr_sinpi (scr.a, scr.a, MPFR_RNDZ);
  round_all (scr.a, inex, mask, out);
}

void
ref_pow (double x, double y, unsigned mask, double out[REF_NRND])
{
  scratch_init ();
  mpfr_set_d (scr.b, x, MPFR_RNDN);
  mpfr_set_d (scr.c, y, MPFR_RNDN);
  int inex = mpfr_pow (scr.a, scr.b, scr.c, MPFR_RNDZ);
  round_all (scr.a, inex, mask, out);
}

void
ref_pown (double x, long long int y, unsigned mask, double out[REF_NRND])
{
  scratch_init ();
  mpfr_set_d (scr.b, x, MPFR_RNDN);
  int inex = mpfr_pown (scr.a, scr.b, y, MPFR_RNDZ);
  round_all (scr.a, inex, mask, out);
}

void
ref_powr (double x, double y, unsigned mask, double out[REF_NRND])
{
  scratch_init ();
  mpfr_set_d (scr.b, x, MPFR_RNDN);
  mpfr_set_d (scr.c, y, MPFR_RNDN);
  int inex = mpfr_powr (scr.a, scr.b, scr.c, MPFR_RNDZ);
  round_all (scr.a, inex, mask, out);
}

void
ref_rootn (double x, long long int y, unsigned mask, double out[REF_NRND])
{
  scratch_init ();
  mpfr_set_d (scr.b, x, MPFR_RNDN);
  int inex = mpfr_rootn_si (scr.a, scr.b, y, MPFR_RNDZ);
  round_all (scr.a, inex, mask, out);
}

void
ref_tan (double x, unsigned mask, double out[REF_NRND])
{
  scratch_init ();
  mpfr_set_d (scr.a, x, MPFR_RNDN);
  int inex = mpfr_tan (scr.a, scr.a, MPFR_RNDZ);
  round_all (scr.a, inex, mask, out);
}

void
ref_tanh (double x, unsigned mask, double out[REF_NRND])
{
  scratch_init ();
  mpfr_set_d (scr.a, x, MPFR_RNDN);
  int inex = mpfr_tanh (scr.a, scr.a, MPFR_RNDZ);
  round_all (scr.a, inex, mask, out);
}

void
ref_tanpi (double x, unsigned mask, double out[REF_NRND])
{
  scratch_init ();
  mpfr_set_d (scr.a, x, MPFR_RNDN);
  int inex = mpfr_tanpi (scr.a, scr.a, MPFR_RNDZ);
  round_all (scr.a, inex, mask, out);
}

void
ref_tgamma (double x, unsigned mask, double out[REF_NRND])
{
  double fx = __builtin_floor (x);
  if (fx == x)
    if (x < 0.0)
      {
	broadcast (__builtin_nanf ("12"), mask, out);
	return;
      }
  scratch_init ();
  mpfr_set_d (scr.a, x, MPFR_RNDN);
  int inex = mpfr_gamma (scr.a, scr.a, MPFR_RNDZ);
  round_all (scr.a, inex, mask, out);
}
