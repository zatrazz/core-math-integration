#include <mpfr.h>
#include <math.h>
#include <stdio.h>
#include <float.h>

static mpfr_rnd_t rnd2[] = { MPFR_RNDN, MPFR_RNDZ, MPFR_RNDU, MPFR_RNDD };

static int rnd = 0; /* default is to nearest */

int ref_fesetround(int rounding_mode)
{
  rnd = rounding_mode;
  return 0;
}

void ref_init(void)
{
  mpfr_set_emin (-148);
  mpfr_set_emax (128);
}

float
ref_rsqrt (float x)
{
  /* mpfr_rec_sqrt differs from IEEE 754-2019: IEEE 754-2019 says that
     rsqrt(-0) should give -Inf, whereas mpfr_rec_sqrt(-0) gives +Inf */
  if (x == 0.0f && 1.0f / x < 0.0f)
    return 1.0f / x;
  mpfr_t y;
  mpfr_init2 (y, 24);
  mpfr_set_flt (y, x, MPFR_RNDN);
  mpfr_rec_sqrt (y, y, rnd2[rnd]);
  /* since |x| < 2^128 for non-zero x, we have 2^-64 < 1/sqrt(x),
     thus no underflow can happen, and there is no need to call
     mpfr_subnormalize */
  float ret = mpfr_get_flt (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

extern float cr_rsqrtf(float x);

static void
check_input (float x)
{
  float r_mpfr = ref_rsqrt (x);
  float r_cr = cr_rsqrtf (x);

  printf ("ref=%a cr=%a\n", r_mpfr, r_cr);
}

int main (int argc, char *argv[])
{
  ref_init ();

  check_input (INFINITY);

  return 0;
}
