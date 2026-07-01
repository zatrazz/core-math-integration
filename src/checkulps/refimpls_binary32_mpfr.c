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
// single high-precision result to float for every rounding mode selected in
// MASK.  See refimpls_binary64_mpfr.c for the details of the round-to-odd
// sticky-bit scheme used by set_sticky()/round_all() and the per-thread
// scratch reuse.

enum
{
  INTERNAL_PRECISION = FLT_MANT_DIG + 2
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
// ternary) to float for every rounding mode selected in MASK.
static void
round_all (mpfr_t hi, int inex, unsigned mask, float out[REF_NRND])
{
  set_sticky (hi, inex);
  for (int i = 0; i < REF_NRND; i++)
    if (mask & (1u << i))
      out[i] = mpfr_get_flt (hi, ref_rnd_modes[i]);
}

// Store the mode-independent value V into every slot selected in MASK.  Used
// for special cases (NaN/Inf/zero) whose result does not depend on rounding.
static void
broadcast (float v, unsigned mask, float out[REF_NRND])
{
  for (int i = 0; i < REF_NRND; i++)
    if (mask & (1u << i))
      out[i] = v;
}

void
ref_acosf (float x, unsigned mask, float out[REF_NRND])
{
  scratch_init ();
  mpfr_set_flt (scr.a, x, MPFR_RNDN);
  int inex = mpfr_acos (scr.a, scr.a, MPFR_RNDZ);
  round_all (scr.a, inex, mask, out);
}

void
ref_acoshf (float x, unsigned mask, float out[REF_NRND])
{
  scratch_init ();
  mpfr_set_flt (scr.a, x, MPFR_RNDN);
  int inex = mpfr_acosh (scr.a, scr.a, MPFR_RNDZ);
  round_all (scr.a, inex, mask, out);
}

void
ref_acospif (float x, unsigned mask, float out[REF_NRND])
{
  scratch_init ();
  mpfr_set_flt (scr.a, x, MPFR_RNDN);
  int inex = mpfr_acospi (scr.a, scr.a, MPFR_RNDZ);
  round_all (scr.a, inex, mask, out);
}

void
ref_asinf (float x, unsigned mask, float out[REF_NRND])
{
  scratch_init ();
  mpfr_set_flt (scr.a, x, MPFR_RNDN);
  int inex = mpfr_asin (scr.a, scr.a, MPFR_RNDZ);
  round_all (scr.a, inex, mask, out);
}

void
ref_asinhf (float x, unsigned mask, float out[REF_NRND])
{
  scratch_init ();
  mpfr_set_flt (scr.a, x, MPFR_RNDN);
  int inex = mpfr_asinh (scr.a, scr.a, MPFR_RNDZ);
  round_all (scr.a, inex, mask, out);
}

void
ref_asinpif (float x, unsigned mask, float out[REF_NRND])
{
  scratch_init ();
  mpfr_set_flt (scr.a, x, MPFR_RNDN);
  int inex = mpfr_asinpi (scr.a, scr.a, MPFR_RNDZ);
  round_all (scr.a, inex, mask, out);
}

void
ref_atanf (float x, unsigned mask, float out[REF_NRND])
{
  scratch_init ();
  mpfr_set_flt (scr.a, x, MPFR_RNDN);
  int inex = mpfr_atan (scr.a, scr.a, MPFR_RNDZ);
  round_all (scr.a, inex, mask, out);
}

void
ref_atan2f (float x, float y, unsigned mask, float out[REF_NRND])
{
  scratch_init ();
  mpfr_set_flt (scr.b, x, MPFR_RNDN);
  mpfr_set_flt (scr.c, y, MPFR_RNDN);
  int inex = mpfr_atan2 (scr.a, scr.b, scr.c, MPFR_RNDZ);
  round_all (scr.a, inex, mask, out);
}

void
ref_atanhf (float x, unsigned mask, float out[REF_NRND])
{
  scratch_init ();
  mpfr_set_flt (scr.a, x, MPFR_RNDN);
  int inex = mpfr_atanh (scr.a, scr.a, MPFR_RNDZ);
  round_all (scr.a, inex, mask, out);
}

void
ref_atanpif (float x, unsigned mask, float out[REF_NRND])
{
  scratch_init ();
  mpfr_set_flt (scr.a, x, MPFR_RNDN);
  int inex = mpfr_atanpi (scr.a, scr.a, MPFR_RNDZ);
  round_all (scr.a, inex, mask, out);
}

void
ref_cbrtf (float x, unsigned mask, float out[REF_NRND])
{
  scratch_init ();
  mpfr_set_flt (scr.a, x, MPFR_RNDN);
  int inex = mpfr_cbrt (scr.a, scr.a, MPFR_RNDZ);
  round_all (scr.a, inex, mask, out);
}

void
ref_compoundnf (float x, long long int y, unsigned mask, float out[REF_NRND])
{
  scratch_init ();
  mpfr_set_flt (scr.b, x, MPFR_RNDN);
  int inex = mpfr_compound_si (scr.a, scr.b, y, MPFR_RNDZ);
  round_all (scr.a, inex, mask, out);
}

void
ref_cosf (float x, unsigned mask, float out[REF_NRND])
{
  scratch_init ();
  mpfr_set_flt (scr.a, x, MPFR_RNDN);
  int inex = mpfr_cos (scr.a, scr.a, MPFR_RNDZ);
  round_all (scr.a, inex, mask, out);
}

void
ref_coshf (float x, unsigned mask, float out[REF_NRND])
{
  scratch_init ();
  mpfr_set_flt (scr.a, x, MPFR_RNDN);
  int inex = mpfr_cosh (scr.a, scr.a, MPFR_RNDZ);
  round_all (scr.a, inex, mask, out);
}

void
ref_cospif (float x, unsigned mask, float out[REF_NRND])
{
  scratch_init ();
  mpfr_set_flt (scr.a, x, MPFR_RNDN);
  int inex = mpfr_cospi (scr.a, scr.a, MPFR_RNDZ);
  round_all (scr.a, inex, mask, out);
}

void
ref_erff (float x, unsigned mask, float out[REF_NRND])
{
  scratch_init ();
  mpfr_set_flt (scr.a, x, MPFR_RNDN);
  int inex = mpfr_erf (scr.a, scr.a, MPFR_RNDZ);
  round_all (scr.a, inex, mask, out);
}

void
ref_erfcf (float x, unsigned mask, float out[REF_NRND])
{
  scratch_init ();
  mpfr_set_flt (scr.a, x, MPFR_RNDN);
  int inex = mpfr_erfc (scr.a, scr.a, MPFR_RNDZ);
  round_all (scr.a, inex, mask, out);
}

void
ref_exp10f (float x, unsigned mask, float out[REF_NRND])
{
  scratch_init ();
  mpfr_set_flt (scr.a, x, MPFR_RNDN);
  int inex = mpfr_exp10 (scr.a, scr.a, MPFR_RNDZ);
  round_all (scr.a, inex, mask, out);
}

void
ref_exp10m1f (float x, unsigned mask, float out[REF_NRND])
{
  scratch_init ();
  mpfr_set_flt (scr.a, x, MPFR_RNDN);
  int inex = mpfr_exp10m1 (scr.a, scr.a, MPFR_RNDZ);
  round_all (scr.a, inex, mask, out);
}

void
ref_exp2f (float x, unsigned mask, float out[REF_NRND])
{
  scratch_init ();
  mpfr_set_flt (scr.a, x, MPFR_RNDN);
  int inex = mpfr_exp2 (scr.a, scr.a, MPFR_RNDZ);
  round_all (scr.a, inex, mask, out);
}

void
ref_exp2m1f (float x, unsigned mask, float out[REF_NRND])
{
  scratch_init ();
  mpfr_set_flt (scr.a, x, MPFR_RNDN);
  int inex = mpfr_exp2m1 (scr.a, scr.a, MPFR_RNDZ);
  round_all (scr.a, inex, mask, out);
}

void
ref_expf (float x, unsigned mask, float out[REF_NRND])
{
  scratch_init ();
  mpfr_set_flt (scr.a, x, MPFR_RNDN);
  int inex = mpfr_exp (scr.a, scr.a, MPFR_RNDZ);
  round_all (scr.a, inex, mask, out);
}

void
ref_expm1f (float x, unsigned mask, float out[REF_NRND])
{
  scratch_init ();
  mpfr_set_flt (scr.a, x, MPFR_RNDN);
  int inex = mpfr_expm1 (scr.a, scr.a, MPFR_RNDZ);
  round_all (scr.a, inex, mask, out);
}

void
ref_hypotf (float x, float y, unsigned mask, float out[REF_NRND])
{
  scratch_init ();
  mpfr_set_flt (scr.b, x, MPFR_RNDN);
  mpfr_set_flt (scr.c, y, MPFR_RNDN);
  int inex = mpfr_hypot (scr.a, scr.b, scr.c, MPFR_RNDZ);
  round_all (scr.a, inex, mask, out);
}

void
ref_lgammaf (float x, unsigned mask, float out[REF_NRND])
{
  int sign;
  scratch_init ();
  mpfr_set_flt (scr.a, x, MPFR_RNDN);
  int inex = mpfr_lgamma (scr.a, &sign, scr.a, MPFR_RNDZ);
  round_all (scr.a, inex, mask, out);
}

void
ref_logf (float x, unsigned mask, float out[REF_NRND])
{
  scratch_init ();
  mpfr_set_flt (scr.a, x, MPFR_RNDN);
  int inex = mpfr_log (scr.a, scr.a, MPFR_RNDZ);
  round_all (scr.a, inex, mask, out);
}

void
ref_log1pf (float x, unsigned mask, float out[REF_NRND])
{
  scratch_init ();
  mpfr_set_flt (scr.a, x, MPFR_RNDN);
  int inex = mpfr_log1p (scr.a, scr.a, MPFR_RNDZ);
  round_all (scr.a, inex, mask, out);
}

void
ref_log2f (float x, unsigned mask, float out[REF_NRND])
{
  scratch_init ();
  mpfr_set_flt (scr.a, x, MPFR_RNDN);
  int inex = mpfr_log2 (scr.a, scr.a, MPFR_RNDZ);
  round_all (scr.a, inex, mask, out);
}

void
ref_log2p1f (float x, unsigned mask, float out[REF_NRND])
{
  scratch_init ();
  mpfr_set_flt (scr.a, x, MPFR_RNDN);
  int inex = mpfr_log2p1 (scr.a, scr.a, MPFR_RNDZ);
  round_all (scr.a, inex, mask, out);
}

void
ref_log10f (float x, unsigned mask, float out[REF_NRND])
{
  scratch_init ();
  mpfr_set_flt (scr.a, x, MPFR_RNDN);
  int inex = mpfr_log10 (scr.a, scr.a, MPFR_RNDZ);
  round_all (scr.a, inex, mask, out);
}

void
ref_log10p1f (float x, unsigned mask, float out[REF_NRND])
{
  scratch_init ();
  mpfr_set_flt (scr.a, x, MPFR_RNDN);
  int inex = mpfr_log10p1 (scr.a, scr.a, MPFR_RNDZ);
  round_all (scr.a, inex, mask, out);
}

void
ref_rsqrtf (float x, unsigned mask, float out[REF_NRND])
{
  /* mpfr_rec_sqrt differs from IEEE 754-2019: IEEE 754-2019 says that
     rsqrt(-0) should give -Inf, whereas mpfr_rec_sqrt(-0) gives +Inf */
  if (x == 0.0f && 1.0f / x < 0.0f)
    {
      broadcast (1.0f / x, mask, out);
      return;
    }
  scratch_init ();
  mpfr_set_flt (scr.a, x, MPFR_RNDN);
  int inex = mpfr_rec_sqrt (scr.a, scr.a, MPFR_RNDZ);
  round_all (scr.a, inex, mask, out);
}

void
ref_sinf (float x, unsigned mask, float out[REF_NRND])
{
  scratch_init ();
  mpfr_set_flt (scr.a, x, MPFR_RNDN);
  int inex = mpfr_sin (scr.a, scr.a, MPFR_RNDZ);
  round_all (scr.a, inex, mask, out);
}

void
ref_sincosf (float x, unsigned mask, float sinp[REF_NRND], float cosp[REF_NRND])
{
  scratch_init ();
  mpfr_set_flt (scr.a, x, MPFR_RNDN);
  int inex_sin = mpfr_sin (scr.a, scr.a, MPFR_RNDZ);
  round_all (scr.a, inex_sin, mask, sinp);

  mpfr_set_flt (scr.b, x, MPFR_RNDN);
  int inex_cos = mpfr_cos (scr.b, scr.b, MPFR_RNDZ);
  round_all (scr.b, inex_cos, mask, cosp);
}

void
ref_sinhf (float x, unsigned mask, float out[REF_NRND])
{
  scratch_init ();
  mpfr_set_flt (scr.a, x, MPFR_RNDN);
  int inex = mpfr_sinh (scr.a, scr.a, MPFR_RNDZ);
  round_all (scr.a, inex, mask, out);
}

void
ref_sinpif (float x, unsigned mask, float out[REF_NRND])
{
  scratch_init ();
  mpfr_set_flt (scr.a, x, MPFR_RNDN);
  int inex = mpfr_sinpi (scr.a, scr.a, MPFR_RNDZ);
  round_all (scr.a, inex, mask, out);
}

void
ref_powf (float x, float y, unsigned mask, float out[REF_NRND])
{
  scratch_init ();
  mpfr_set_flt (scr.b, x, MPFR_RNDN);
  mpfr_set_flt (scr.c, y, MPFR_RNDN);
  int inex = mpfr_pow (scr.a, scr.b, scr.c, MPFR_RNDZ);
  round_all (scr.a, inex, mask, out);
}

void
ref_pownf (float x, long long int y, unsigned mask, float out[REF_NRND])
{
  scratch_init ();
  mpfr_set_flt (scr.b, x, MPFR_RNDN);
  int inex = mpfr_pown (scr.a, scr.b, y, MPFR_RNDZ);
  round_all (scr.a, inex, mask, out);
}

void
ref_powrf (float x, float y, unsigned mask, float out[REF_NRND])
{
  scratch_init ();
  mpfr_set_flt (scr.b, x, MPFR_RNDN);
  mpfr_set_flt (scr.c, y, MPFR_RNDN);
  int inex = mpfr_powr (scr.a, scr.b, scr.c, MPFR_RNDZ);
  round_all (scr.a, inex, mask, out);
}

void
ref_rootnf (float x, long long int y, unsigned mask, float out[REF_NRND])
{
  scratch_init ();
  mpfr_set_flt (scr.b, x, MPFR_RNDN);
  int inex = mpfr_rootn_si (scr.a, scr.b, y, MPFR_RNDZ);
  round_all (scr.a, inex, mask, out);
}

void
ref_tanf (float x, unsigned mask, float out[REF_NRND])
{
  scratch_init ();
  mpfr_set_flt (scr.a, x, MPFR_RNDN);
  int inex = mpfr_tan (scr.a, scr.a, MPFR_RNDZ);
  round_all (scr.a, inex, mask, out);
}

void
ref_tanhf (float x, unsigned mask, float out[REF_NRND])
{
  scratch_init ();
  mpfr_set_flt (scr.a, x, MPFR_RNDN);
  int inex = mpfr_tanh (scr.a, scr.a, MPFR_RNDZ);
  round_all (scr.a, inex, mask, out);
}

void
ref_tanpif (float x, unsigned mask, float out[REF_NRND])
{
  scratch_init ();
  mpfr_set_flt (scr.a, x, MPFR_RNDN);
  int inex = mpfr_tanpi (scr.a, scr.a, MPFR_RNDZ);
  round_all (scr.a, inex, mask, out);
}

void
ref_tgammaf (float x, unsigned mask, float out[REF_NRND])
{
  scratch_init ();
  mpfr_set_flt (scr.a, x, MPFR_RNDN);
  int inex = mpfr_gamma (scr.a, scr.a, MPFR_RNDZ);
  round_all (scr.a, inex, mask, out);
}
