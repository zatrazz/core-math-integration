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
// sticky-bit scheme used by set_sticky()/round_all().

enum
{
  INTERNAL_PRECISION = FLT_MANT_DIG + 2
};

static const mpfr_rnd_t ref_rnd_modes[REF_NRND]
    = { MPFR_RNDN, MPFR_RNDU, MPFR_RNDD, MPFR_RNDZ };

// Fold the round-toward-zero ternary INEX into a sticky bit (round-to-odd):
// if the result was inexact, force the least significant bit of the mantissa.
static void
set_sticky (mpfr_t r, int inex)
{
  if (inex == 0 || !mpfr_number_p (r))
    return;
  mpz_t m;
  mpz_init (m);
  mpfr_exp_t e = mpfr_get_z_2exp (m, r);
  if (mpz_sgn (m) < 0)
    {
      mpz_neg (m, m);
      mpz_setbit (m, 0);
      mpz_neg (m, m);
    }
  else
    mpz_setbit (m, 0);
  mpfr_set_z_2exp (r, m, e, MPFR_RNDN);
  mpz_clear (m);
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
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_flt (y, x, MPFR_RNDN);
  int inex = mpfr_acos (y, y, MPFR_RNDZ);
  round_all (y, inex, mask, out);
  mpfr_clear (y);
}

void
ref_acoshf (float x, unsigned mask, float out[REF_NRND])
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_flt (y, x, MPFR_RNDN);
  int inex = mpfr_acosh (y, y, MPFR_RNDZ);
  round_all (y, inex, mask, out);
  mpfr_clear (y);
}

void
ref_acospif (float x, unsigned mask, float out[REF_NRND])
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_flt (y, x, MPFR_RNDN);
  int inex = mpfr_acospi (y, y, MPFR_RNDZ);
  round_all (y, inex, mask, out);
  mpfr_clear (y);
}

void
ref_asinf (float x, unsigned mask, float out[REF_NRND])
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_flt (y, x, MPFR_RNDN);
  int inex = mpfr_asin (y, y, MPFR_RNDZ);
  round_all (y, inex, mask, out);
  mpfr_clear (y);
}

void
ref_asinhf (float x, unsigned mask, float out[REF_NRND])
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_flt (y, x, MPFR_RNDN);
  int inex = mpfr_asinh (y, y, MPFR_RNDZ);
  round_all (y, inex, mask, out);
  mpfr_clear (y);
}

void
ref_asinpif (float x, unsigned mask, float out[REF_NRND])
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_flt (y, x, MPFR_RNDN);
  int inex = mpfr_asinpi (y, y, MPFR_RNDZ);
  round_all (y, inex, mask, out);
  mpfr_clear (y);
}

void
ref_atanf (float x, unsigned mask, float out[REF_NRND])
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_flt (y, x, MPFR_RNDN);
  int inex = mpfr_atan (y, y, MPFR_RNDZ);
  round_all (y, inex, mask, out);
  mpfr_clear (y);
}

void
ref_atan2f (float x, float y, unsigned mask, float out[REF_NRND])
{
  mpfr_t xm, ym, zm;
  mpfr_init2 (xm, INTERNAL_PRECISION);
  mpfr_init2 (ym, INTERNAL_PRECISION);
  mpfr_init2 (zm, INTERNAL_PRECISION);
  mpfr_set_flt (xm, x, MPFR_RNDN);
  mpfr_set_flt (ym, y, MPFR_RNDN);
  int inex = mpfr_atan2 (zm, xm, ym, MPFR_RNDZ);
  round_all (zm, inex, mask, out);
  mpfr_clear (xm);
  mpfr_clear (ym);
  mpfr_clear (zm);
}

void
ref_atanhf (float x, unsigned mask, float out[REF_NRND])
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_flt (y, x, MPFR_RNDN);
  int inex = mpfr_atanh (y, y, MPFR_RNDZ);
  round_all (y, inex, mask, out);
  mpfr_clear (y);
}

void
ref_atanpif (float x, unsigned mask, float out[REF_NRND])
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_flt (y, x, MPFR_RNDN);
  int inex = mpfr_atanpi (y, y, MPFR_RNDZ);
  round_all (y, inex, mask, out);
  mpfr_clear (y);
}

void
ref_cbrtf (float x, unsigned mask, float out[REF_NRND])
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_flt (y, x, MPFR_RNDN);
  int inex = mpfr_cbrt (y, y, MPFR_RNDZ);
  round_all (y, inex, mask, out);
  mpfr_clear (y);
}

void
ref_compoundnf (float x, long long int y, unsigned mask, float out[REF_NRND])
{
  mpfr_t xm, zm;
  mpfr_init2 (xm, INTERNAL_PRECISION);
  mpfr_init2 (zm, INTERNAL_PRECISION);
  mpfr_set_flt (xm, x, MPFR_RNDN);
  int inex = mpfr_compound_si (zm, xm, y, MPFR_RNDZ);
  round_all (zm, inex, mask, out);
  mpfr_clear (xm);
  mpfr_clear (zm);
}

void
ref_cosf (float x, unsigned mask, float out[REF_NRND])
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_flt (y, x, MPFR_RNDN);
  int inex = mpfr_cos (y, y, MPFR_RNDZ);
  round_all (y, inex, mask, out);
  mpfr_clear (y);
}

void
ref_coshf (float x, unsigned mask, float out[REF_NRND])
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_flt (y, x, MPFR_RNDN);
  int inex = mpfr_cosh (y, y, MPFR_RNDZ);
  round_all (y, inex, mask, out);
  mpfr_clear (y);
}

void
ref_cospif (float x, unsigned mask, float out[REF_NRND])
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_flt (y, x, MPFR_RNDN);
  int inex = mpfr_cospi (y, y, MPFR_RNDZ);
  round_all (y, inex, mask, out);
  mpfr_clear (y);
}

void
ref_erff (float x, unsigned mask, float out[REF_NRND])
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_flt (y, x, MPFR_RNDN);
  int inex = mpfr_erf (y, y, MPFR_RNDZ);
  round_all (y, inex, mask, out);
  mpfr_clear (y);
}

void
ref_erfcf (float x, unsigned mask, float out[REF_NRND])
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_flt (y, x, MPFR_RNDN);
  int inex = mpfr_erfc (y, y, MPFR_RNDZ);
  round_all (y, inex, mask, out);
  mpfr_clear (y);
}

void
ref_exp10f (float x, unsigned mask, float out[REF_NRND])
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_flt (y, x, MPFR_RNDN);
  int inex = mpfr_exp10 (y, y, MPFR_RNDZ);
  round_all (y, inex, mask, out);
  mpfr_clear (y);
}

void
ref_exp10m1f (float x, unsigned mask, float out[REF_NRND])
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_flt (y, x, MPFR_RNDN);
  int inex = mpfr_exp10m1 (y, y, MPFR_RNDZ);
  round_all (y, inex, mask, out);
  mpfr_clear (y);
}

void
ref_exp2f (float x, unsigned mask, float out[REF_NRND])
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_flt (y, x, MPFR_RNDN);
  int inex = mpfr_exp2 (y, y, MPFR_RNDZ);
  round_all (y, inex, mask, out);
  mpfr_clear (y);
}

void
ref_exp2m1f (float x, unsigned mask, float out[REF_NRND])
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_flt (y, x, MPFR_RNDN);
  int inex = mpfr_exp2m1 (y, y, MPFR_RNDZ);
  round_all (y, inex, mask, out);
  mpfr_clear (y);
}

void
ref_expf (float x, unsigned mask, float out[REF_NRND])
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_flt (y, x, MPFR_RNDN);
  int inex = mpfr_exp (y, y, MPFR_RNDZ);
  round_all (y, inex, mask, out);
  mpfr_clear (y);
}

void
ref_expm1f (float x, unsigned mask, float out[REF_NRND])
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_flt (y, x, MPFR_RNDN);
  int inex = mpfr_expm1 (y, y, MPFR_RNDZ);
  round_all (y, inex, mask, out);
  mpfr_clear (y);
}

void
ref_hypotf (float x, float y, unsigned mask, float out[REF_NRND])
{
  mpfr_t xm, ym, zm;
  mpfr_init2 (xm, INTERNAL_PRECISION);
  mpfr_init2 (ym, INTERNAL_PRECISION);
  mpfr_init2 (zm, INTERNAL_PRECISION);
  mpfr_set_flt (xm, x, MPFR_RNDN);
  mpfr_set_flt (ym, y, MPFR_RNDN);
  int inex = mpfr_hypot (zm, xm, ym, MPFR_RNDZ);
  round_all (zm, inex, mask, out);
  mpfr_clear (xm);
  mpfr_clear (ym);
  mpfr_clear (zm);
}

void
ref_lgammaf (float x, unsigned mask, float out[REF_NRND])
{
  int sign;
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_flt (y, x, MPFR_RNDN);
  int inex = mpfr_lgamma (y, &sign, y, MPFR_RNDZ);
  round_all (y, inex, mask, out);
  mpfr_clear (y);
}

void
ref_logf (float x, unsigned mask, float out[REF_NRND])
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_flt (y, x, MPFR_RNDN);
  int inex = mpfr_log (y, y, MPFR_RNDZ);
  round_all (y, inex, mask, out);
  mpfr_clear (y);
}

void
ref_log1pf (float x, unsigned mask, float out[REF_NRND])
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_flt (y, x, MPFR_RNDN);
  int inex = mpfr_log1p (y, y, MPFR_RNDZ);
  round_all (y, inex, mask, out);
  mpfr_clear (y);
}

void
ref_log2f (float x, unsigned mask, float out[REF_NRND])
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_flt (y, x, MPFR_RNDN);
  int inex = mpfr_log2 (y, y, MPFR_RNDZ);
  round_all (y, inex, mask, out);
  mpfr_clear (y);
}

void
ref_log2p1f (float x, unsigned mask, float out[REF_NRND])
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_flt (y, x, MPFR_RNDN);
  int inex = mpfr_log2p1 (y, y, MPFR_RNDZ);
  round_all (y, inex, mask, out);
  mpfr_clear (y);
}

void
ref_log10f (float x, unsigned mask, float out[REF_NRND])
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_flt (y, x, MPFR_RNDN);
  int inex = mpfr_log10 (y, y, MPFR_RNDZ);
  round_all (y, inex, mask, out);
  mpfr_clear (y);
}

void
ref_log10p1f (float x, unsigned mask, float out[REF_NRND])
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_flt (y, x, MPFR_RNDN);
  int inex = mpfr_log10p1 (y, y, MPFR_RNDZ);
  round_all (y, inex, mask, out);
  mpfr_clear (y);
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
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_flt (y, x, MPFR_RNDN);
  int inex = mpfr_rec_sqrt (y, y, MPFR_RNDZ);
  round_all (y, inex, mask, out);
  mpfr_clear (y);
}

void
ref_sinf (float x, unsigned mask, float out[REF_NRND])
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_flt (y, x, MPFR_RNDN);
  int inex = mpfr_sin (y, y, MPFR_RNDZ);
  round_all (y, inex, mask, out);
  mpfr_clear (y);
}

void
ref_sincosf (float x, unsigned mask, float sinp[REF_NRND], float cosp[REF_NRND])
{
  mpfr_t xm0, xm1;
  mpfr_init2 (xm0, INTERNAL_PRECISION);
  mpfr_init2 (xm1, INTERNAL_PRECISION);

  mpfr_set_flt (xm0, x, MPFR_RNDN);
  int inex_sin = mpfr_sin (xm0, xm0, MPFR_RNDZ);
  round_all (xm0, inex_sin, mask, sinp);

  mpfr_set_flt (xm1, x, MPFR_RNDN);
  int inex_cos = mpfr_cos (xm1, xm1, MPFR_RNDZ);
  round_all (xm1, inex_cos, mask, cosp);

  mpfr_clear (xm0);
  mpfr_clear (xm1);
}

void
ref_sinhf (float x, unsigned mask, float out[REF_NRND])
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_flt (y, x, MPFR_RNDN);
  int inex = mpfr_sinh (y, y, MPFR_RNDZ);
  round_all (y, inex, mask, out);
  mpfr_clear (y);
}

void
ref_sinpif (float x, unsigned mask, float out[REF_NRND])
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_flt (y, x, MPFR_RNDN);
  int inex = mpfr_sinpi (y, y, MPFR_RNDZ);
  round_all (y, inex, mask, out);
  mpfr_clear (y);
}

void
ref_powf (float x, float y, unsigned mask, float out[REF_NRND])
{
  mpfr_t xm, ym, zm;
  mpfr_init2 (xm, INTERNAL_PRECISION);
  mpfr_init2 (ym, INTERNAL_PRECISION);
  mpfr_init2 (zm, INTERNAL_PRECISION);
  mpfr_set_flt (xm, x, MPFR_RNDN);
  mpfr_set_flt (ym, y, MPFR_RNDN);
  int inex = mpfr_pow (zm, xm, ym, MPFR_RNDZ);
  round_all (zm, inex, mask, out);
  mpfr_clear (xm);
  mpfr_clear (ym);
  mpfr_clear (zm);
}

void
ref_pownf (float x, long long int y, unsigned mask, float out[REF_NRND])
{
  mpfr_t xm, zm;
  mpfr_init2 (xm, INTERNAL_PRECISION);
  mpfr_init2 (zm, INTERNAL_PRECISION);
  mpfr_set_flt (xm, x, MPFR_RNDN);
  int inex = mpfr_pown (zm, xm, y, MPFR_RNDZ);
  round_all (zm, inex, mask, out);
  mpfr_clear (xm);
  mpfr_clear (zm);
}

void
ref_powrf (float x, float y, unsigned mask, float out[REF_NRND])
{
  mpfr_t xm, ym, zm;
  mpfr_init2 (xm, INTERNAL_PRECISION);
  mpfr_init2 (ym, INTERNAL_PRECISION);
  mpfr_init2 (zm, INTERNAL_PRECISION);
  mpfr_set_flt (xm, x, MPFR_RNDN);
  mpfr_set_flt (ym, y, MPFR_RNDN);
  int inex = mpfr_powr (zm, xm, ym, MPFR_RNDZ);
  round_all (zm, inex, mask, out);
  mpfr_clear (xm);
  mpfr_clear (ym);
  mpfr_clear (zm);
}

void
ref_rootnf (float x, long long int y, unsigned mask, float out[REF_NRND])
{
  mpfr_t xm, zm;
  mpfr_init2 (xm, INTERNAL_PRECISION);
  mpfr_init2 (zm, INTERNAL_PRECISION);
  mpfr_set_flt (xm, x, MPFR_RNDN);
  int inex = mpfr_rootn_si (zm, xm, y, MPFR_RNDZ);
  round_all (zm, inex, mask, out);
  mpfr_clear (xm);
  mpfr_clear (zm);
}

void
ref_tanf (float x, unsigned mask, float out[REF_NRND])
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_flt (y, x, MPFR_RNDN);
  int inex = mpfr_tan (y, y, MPFR_RNDZ);
  round_all (y, inex, mask, out);
  mpfr_clear (y);
}

void
ref_tanhf (float x, unsigned mask, float out[REF_NRND])
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_flt (y, x, MPFR_RNDN);
  int inex = mpfr_tanh (y, y, MPFR_RNDZ);
  round_all (y, inex, mask, out);
  mpfr_clear (y);
}

void
ref_tanpif (float x, unsigned mask, float out[REF_NRND])
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_flt (y, x, MPFR_RNDN);
  int inex = mpfr_tanpi (y, y, MPFR_RNDZ);
  round_all (y, inex, mask, out);
  mpfr_clear (y);
}

void
ref_tgammaf (float x, unsigned mask, float out[REF_NRND])
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_flt (y, x, MPFR_RNDN);
  int inex = mpfr_gamma (y, y, MPFR_RNDZ);
  round_all (y, inex, mask, out);
  mpfr_clear (y);
}
