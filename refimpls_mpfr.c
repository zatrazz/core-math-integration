#include <stdint.h>
#include <mpfr.h>
#include <math.h>

double
ref_acos (double x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_init2 (y, 53);
  mpfr_set_d (y, x, MPFR_RNDN);
  int inex = mpfr_acos (y, y, rnd);
  mpfr_subnormalize (y, inex, rnd);
  double ret = mpfr_get_d (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

/* code from MPFR */
double
ref_acospi (double x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_init2 (y, 53);
  mpfr_set_d (y, x, MPFR_RNDN);
  mpfr_acospi (y, y, rnd);
  /* no need to call mpfr_subnormalize since the smallest non-zero value
     is obtained for x=0x1.fffffffffffffp-1, and is 0x1.45f306dc9c883p-28
     (rounded to nearest) */
  double ret = mpfr_get_d (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

typedef union {double f; uint64_t u;} b64u64_u;

double
ref_acosh (double x, mpfr_rnd_t rnd)
{
  b64u64_u ix = {.f = x};
  if(__builtin_expect(ix.u<=0x3ff0000000000000ull, 0)){
    if(ix.u==0x3ff0000000000000ull) return 0;
    return __builtin_nan("x<1");
  }
  if(__builtin_expect(ix.u>=0x7ff0000000000000ull, 0)){
    uint64_t aix = ix.u<<1;
    if(ix.u==0x7ff0000000000000ull || aix>((uint64_t)0x7ff<<53)) return x; // +inf or nan
    return __builtin_nan("x<1");
  }

  mpfr_t y;
  mpfr_init2 (y, 53);
  mpfr_set_d (y, x, MPFR_RNDN);
  mpfr_acosh (y, y, rnd);
  double ret = mpfr_get_d(y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

double
ref_asin (double x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_init2 (y, 53);
  mpfr_set_d (y, x, MPFR_RNDN);
  int inex = mpfr_asin (y, y, rnd);
  mpfr_subnormalize (y, inex, rnd);
  double ret = mpfr_get_d (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

double
ref_asinpi (double x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_init2 (y, 53);
  mpfr_set_d (y, x, MPFR_RNDN);
  int inex = mpfr_asinpi (y, y, rnd);
  mpfr_subnormalize (y, inex, rnd);
  double ret = mpfr_get_d (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

double ref_asinh (double x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_init2(y, 53);
  mpfr_set_d(y, x, MPFR_RNDN);
  int d = mpfr_asinh(y, y, rnd);
  mpfr_subnormalize(y, d, rnd);
  double ret = mpfr_get_d(y, MPFR_RNDN);
  mpfr_clear(y);
  return ret;
}

double ref_atan (double x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_init2 (y, 53);
  mpfr_set_d (y, x, MPFR_RNDN);
  int inex = mpfr_atan (y, y, rnd);
  mpfr_subnormalize (y, inex, rnd);
  double ret = mpfr_get_d (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

double ref_atan2 (double y, double x, mpfr_rnd_t rnd)
{
  mpfr_t z, _x, _y;
  mpfr_inits2(53, z, _x, _y, NULL);
  mpfr_set_d(_x, x, MPFR_RNDN);
  mpfr_set_d(_y, y, MPFR_RNDN);
  int inex = mpfr_atan2 (z, _y, _x, rnd);
  mpfr_subnormalize(z, inex, rnd);
  double ret = mpfr_get_d(z, rnd);
  mpfr_clears(z, _x, _y, NULL);
  return ret;
}

double ref_atanh (double x,  mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_init2(y, 53);
  mpfr_set_d(y, x, MPFR_RNDN);
  int inex = mpfr_atanh(y, y, rnd);
  mpfr_subnormalize(y, inex, rnd);
  double r = mpfr_get_d(y, MPFR_RNDN);
  mpfr_clear(y);
  return r;
}

double
ref_atanpi (double x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_init2 (y, 53);
  mpfr_set_d (y, x, MPFR_RNDN);
  int inex = mpfr_atanpi (y, y, rnd);
  mpfr_subnormalize (y, inex, rnd);
  double ret = mpfr_get_d (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

double ref_sin (double x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_init2 (y, 53);
  mpfr_set_d (y, x, MPFR_RNDN);
  int inex = mpfr_sin (y, y, rnd);
  mpfr_subnormalize (y, inex, rnd);
  double ret = mpfr_get_d (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

double ref_cos (double x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_init2 (y, 53);
  mpfr_set_d (y, x, MPFR_RNDN);
  int inex = mpfr_cos (y, y, rnd);
  mpfr_subnormalize (y, inex, rnd);
  double ret = mpfr_get_d (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

double ref_tan (double x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_init2 (y, 53);
  mpfr_set_d (y, x, MPFR_RNDN);
  int inex = mpfr_tan (y, y, rnd);
  mpfr_subnormalize (y, inex, rnd);
  double ret = mpfr_get_d (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

double
ref_tanh (double x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_init2 (y, 53);
  mpfr_set_d (y, x, MPFR_RNDN);
  int inex = mpfr_tanh (y, y, rnd);
  mpfr_subnormalize (y, inex, rnd);
  double ret = mpfr_get_d (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

double ref_tanpi (double x, mpfr_rnd_t rnd)
{
  if(isnan(x)) return x;
  if(isinf(x)) return -__builtin_nan("");
  mpfr_t y;
  mpfr_init2(y, 53);
  mpfr_set_d(y, x, MPFR_RNDN);
  int inex = mpfr_tanpi(y, y, rnd);
  mpfr_subnormalize(y, inex, rnd);
  double ret = mpfr_get_d(y, MPFR_RNDN);
  mpfr_clear(y);
  return ret;
}
