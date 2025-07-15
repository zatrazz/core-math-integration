#include <stdint.h>
#include <mpfr.h>
#include <math.h>

float
ref_atanf (float x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_init2 (y, 24);
  mpfr_set_flt (y, x, MPFR_RNDN);
  int inex = mpfr_atan (y, y, rnd);
  mpfr_subnormalize (y, inex, rnd);
  float ret = mpfr_get_flt (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

float
ref_atanpif (float x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_init2 (y, 24);
  mpfr_set_flt (y, x, MPFR_RNDN);
  int inex = mpfr_atanpi (y, y, rnd);
  mpfr_subnormalize (y, inex, rnd);
  float ret = mpfr_get_flt (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

float
ref_powf (float x, float y, mpfr_rnd_t rnd)
{
  mpfr_t xm, ym, zm;
  mpfr_init2 (xm, 24);
  mpfr_init2 (ym, 24);
  mpfr_init2 (zm, 24);
  mpfr_set_flt (xm, x, MPFR_RNDN);
  mpfr_set_flt (ym, y, MPFR_RNDN);
  int inex = mpfr_pow (zm, xm, ym, rnd);
  mpfr_subnormalize (zm, inex, rnd);
  float ret = mpfr_get_flt (zm, MPFR_RNDN);
  mpfr_clear (xm);
  mpfr_clear (ym);
  mpfr_clear (zm);
  return ret;
}
