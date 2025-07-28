#include <float.h>
#include <math.h>
#include <mpfr.h>
#include <stdint.h>

enum
{
  INTERNAL_PRECISION = DBL_MANT_DIG
};

double
ref_acos (double x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_d (y, x, MPFR_RNDN);
  int inex = mpfr_acos (y, y, rnd);
  mpfr_subnormalize (y, inex, rnd);
  double ret = mpfr_get_d (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

double
ref_acosh (double x, mpfr_rnd_t rnd)
{
  union
  {
    uint64_t u;
    double f;
  } ix = { .f = x };
  if (__builtin_expect (ix.u <= 0x3ff0000000000000ull, 0))
    {
      if (ix.u == 0x3ff0000000000000ull)
        return 0;
      return __builtin_nan ("x<1");
    }
  if (__builtin_expect (ix.u >= 0x7ff0000000000000ull, 0))
    {
      uint64_t aix = ix.u << 1;
      if (ix.u == 0x7ff0000000000000ull
          || aix > ((uint64_t)0x7ff << INTERNAL_PRECISION))
        return x; // +inf or nan
      return __builtin_nan ("x<1");
    }

  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_d (y, x, MPFR_RNDN);
  mpfr_acosh (y, y, rnd);
  double ret = mpfr_get_d (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

double
ref_acospi (double x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_d (y, x, MPFR_RNDN);
  mpfr_acospi (y, y, rnd);
  /* no need to call mpfr_subnormalize since the smallest non-zero value
     is obtained for x=0x1.fffffffffffffp-1, and is 0x1.45f306dc9c883p-28
     (rounded to nearest) */
  double ret = mpfr_get_d (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

double
ref_asin (double x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_d (y, x, MPFR_RNDN);
  int inex = mpfr_asin (y, y, rnd);
  mpfr_subnormalize (y, inex, rnd);
  double ret = mpfr_get_d (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

double
ref_asinh (double x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_d (y, x, MPFR_RNDN);
  int d = mpfr_asinh (y, y, rnd);
  mpfr_subnormalize (y, d, rnd);
  double ret = mpfr_get_d (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

double
ref_asinpi (double x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_d (y, x, MPFR_RNDN);
  int inex = mpfr_asinpi (y, y, rnd);
  mpfr_subnormalize (y, inex, rnd);
  double ret = mpfr_get_d (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

double
ref_atan (double x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_d (y, x, MPFR_RNDN);
  int inex = mpfr_atan (y, y, rnd);
  mpfr_subnormalize (y, inex, rnd);
  double ret = mpfr_get_d (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

double
ref_atan2 (double y, double x, mpfr_rnd_t rnd)
{
  mpfr_t z, _x, _y;
  mpfr_inits2 (INTERNAL_PRECISION, z, _x, _y, NULL);
  mpfr_set_d (_x, x, MPFR_RNDN);
  mpfr_set_d (_y, y, MPFR_RNDN);
  int inex = mpfr_atan2 (z, _y, _x, rnd);
  mpfr_subnormalize (z, inex, rnd);
  double ret = mpfr_get_d (z, rnd);
  mpfr_clears (z, _x, _y, NULL);
  return ret;
}

double
ref_atanh (double x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_d (y, x, MPFR_RNDN);
  int inex = mpfr_atanh (y, y, rnd);
  mpfr_subnormalize (y, inex, rnd);
  double r = mpfr_get_d (y, MPFR_RNDN);
  mpfr_clear (y);
  return r;
}

double
ref_atanpi (double x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_d (y, x, MPFR_RNDN);
  int inex = mpfr_atanpi (y, y, rnd);
  mpfr_subnormalize (y, inex, rnd);
  double ret = mpfr_get_d (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

double
ref_cbrt (double x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_d (y, x, MPFR_RNDN);
  int inex = mpfr_cbrt (y, y, rnd);
  mpfr_subnormalize (y, inex, rnd);
  double ret = mpfr_get_d (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

double
ref_compoundn (double x, long long int y, mpfr_rnd_t rnd)
{
  mpfr_t xm, zm;
  mpfr_init2 (xm, INTERNAL_PRECISION);
  mpfr_init2 (zm, INTERNAL_PRECISION);
  mpfr_set_d (xm, x, MPFR_RNDN);
  int inex = mpfr_compound_si (zm, xm, y, rnd);
  mpfr_subnormalize (zm, inex, rnd);
  double ret = mpfr_get_d (zm, MPFR_RNDN);
  mpfr_clear (xm);
  mpfr_clear (zm);
  return ret;
}

double
ref_cos (double x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_d (y, x, MPFR_RNDN);
  int inex = mpfr_cos (y, y, rnd);
  mpfr_subnormalize (y, inex, rnd);
  double ret = mpfr_get_d (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

double
ref_cosh (double x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_d (y, x, MPFR_RNDN);
  mpfr_cosh (y, y, rnd);
  /* no need to call mpfr_subnormalize(), since cosh(x) >= 1 */
  double ret = mpfr_get_d (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

double
ref_cospi (double x, mpfr_rnd_t rnd)
{
  if (isnan (x))
    return x;
  if (isinf (x))
    return __builtin_nan ("");
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_d (y, x, MPFR_RNDN);
  mpfr_cospi (y, y, rnd);
  double r = mpfr_get_d (y, MPFR_RNDN);
  mpfr_clear (y);
  return r;
}

double
ref_erf (double x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_d (y, x, MPFR_RNDN);
  int inex = mpfr_erf (y, y, rnd);
  mpfr_subnormalize (y, inex, rnd);
  double ret = mpfr_get_d (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

double
ref_erfc (double x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_d (y, x, MPFR_RNDN);
  int inex = mpfr_erfc (y, y, rnd);
  mpfr_subnormalize (y, inex, rnd);
  double ret = mpfr_get_d (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

double
ref_exp10 (double x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_exp_t emin = mpfr_get_emin ();
  mpfr_set_emin (-1073);
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_d (y, x, MPFR_RNDN);
  int inex = mpfr_exp10 (y, y, rnd);
  mpfr_subnormalize (y, inex, rnd);
  double ret = mpfr_get_d (y, MPFR_RNDN);
  mpfr_clear (y);
  mpfr_set_emin (emin);
  return ret;
}

double
ref_exp10m1 (double x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_d (y, x, MPFR_RNDN);
  int inex = mpfr_exp10m1 (y, y, rnd);
  mpfr_subnormalize (y, inex, rnd);
  double ret = mpfr_get_d (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

double
ref_exp2 (double x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_exp_t emin = mpfr_get_emin ();
  mpfr_set_emin (-1073);
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_d (y, x, MPFR_RNDN);
  int inex = mpfr_exp2 (y, y, rnd);
  mpfr_subnormalize (y, inex, rnd);
  double ret = mpfr_get_d (y, MPFR_RNDN);
  mpfr_clear (y);
  mpfr_set_emin (emin);
  return ret;
}

double
ref_exp2m1 (double x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_d (y, x, MPFR_RNDN);
  int inex = mpfr_exp2m1 (y, y, rnd);
  mpfr_subnormalize (y, inex, rnd);
  double ret = mpfr_get_d (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

double
ref_exp (double x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_exp_t emin = mpfr_get_emin ();
  mpfr_set_emin (-1073);
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_d (y, x, MPFR_RNDN);
  int inex = mpfr_exp (y, y, rnd);
  mpfr_subnormalize (y, inex, rnd);
  double ret = mpfr_get_d (y, MPFR_RNDN);
  mpfr_clear (y);
  mpfr_set_emin (emin);
  return ret;
}

double
ref_expm1 (double x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_d (y, x, MPFR_RNDN);
  int inex = mpfr_expm1 (y, y, rnd);
  mpfr_subnormalize (y, inex, rnd);
  double ret = mpfr_get_d (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

double
ref_hypot (double x, double y, mpfr_rnd_t rnd)
{
  /* since MPFR does not distinguish between quiet/signaling NaN,
     we have to deal with them separately to apply the IEEE rules */
  union
  {
    double f;
    uint64_t u;
  } xi = { .f = x }, yi = { .f = y };
  if ((xi.u << 1) < (0xfffull << 52)
      && (xi.u << 1) > (0x7ffull << INTERNAL_PRECISION)) // x = sNAN
    return x + y;                                        // will return qNAN
  if ((yi.u << 1) < (0xfffull << 52)
      && (yi.u << 1) > (0x7ffull << INTERNAL_PRECISION)) // y = sNAN
    return x + y;                                        // will return qNAN
  if ((xi.u << 1) == 0)
    { // x = +/-0
      yi.u = (yi.u << 1) >> 1;
      return yi.f;
    }
  if ((yi.u << 1) == 0)
    { // y = +/-0
      xi.u = (xi.u << 1) >> 1;
      return xi.f;
    }

  mpfr_t xm, ym, zm;
  mpfr_set_emin (-1073);
  mpfr_init2 (xm, INTERNAL_PRECISION);
  mpfr_init2 (ym, INTERNAL_PRECISION);
  mpfr_init2 (zm, INTERNAL_PRECISION);
  mpfr_set_d (xm, x, MPFR_RNDN);
  mpfr_set_d (ym, y, MPFR_RNDN);
  int inex = mpfr_hypot (zm, xm, ym, rnd);
  mpfr_subnormalize (zm, inex, rnd);
  double ret = mpfr_get_d (zm, MPFR_RNDN);
  mpfr_clear (xm);
  mpfr_clear (ym);
  mpfr_clear (zm);
  return ret;
}

double
ref_lgamma (double x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_d (y, x, MPFR_RNDN);
  int inex = mpfr_lgamma (y, &signgam, y, rnd);
  mpfr_subnormalize (y, inex, rnd);
  double ret = mpfr_get_d (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

double
ref_log (double x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_d (y, x, MPFR_RNDN);
  mpfr_log (y, y, rnd);
  double ret = mpfr_get_d (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

double
ref_log1p (double x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_d (y, x, MPFR_RNDN);
  int inex = mpfr_log1p (y, y, rnd);
  mpfr_subnormalize (y, inex, rnd);
  double ret = mpfr_get_d (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

double
ref_log2 (double x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_d (y, x, MPFR_RNDN);
  mpfr_log2 (y, y, rnd);
  double ret = mpfr_get_d (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

double
ref_log2p1 (double x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_d (y, x, MPFR_RNDN);
  int inex = mpfr_log2p1 (y, y, rnd);
  mpfr_subnormalize (y, inex, rnd);
  double ret = mpfr_get_d (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

double
ref_log10 (double x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_d (y, x, MPFR_RNDN);
  mpfr_log10 (y, y, rnd);
  double ret = mpfr_get_d (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

double
ref_log10p1 (double x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_d (y, x, MPFR_RNDN);
  int inex = mpfr_log10p1 (y, y, rnd);
  mpfr_subnormalize (y, inex, rnd);
  double ret = mpfr_get_d (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

double
ref_rsqrt (double x, mpfr_rnd_t rnd)
{
  /* mpfr_rec_sqrt differs from IEEE 754-2019: IEEE 754-2019 says that
     rsqrt(-0) should give -Inf, whereas mpfr_rec_sqrt(-0) gives +Inf */
  if (x == 0.0 && 1.0 / x < 0)
    return 1.0 / x;
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_d (y, x, MPFR_RNDN);
  mpfr_rec_sqrt (y, y, rnd);
  double ret = mpfr_get_d (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

double
ref_sin (double x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_d (y, x, MPFR_RNDN);
  int inex = mpfr_sin (y, y, rnd);
  mpfr_subnormalize (y, inex, rnd);
  double ret = mpfr_get_d (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

double
ref_sinh (double x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_d (y, x, MPFR_RNDN);
  int inex = mpfr_sinh (y, y, rnd);
  mpfr_subnormalize (y, inex, rnd);
  double ret = mpfr_get_d (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

double
ref_sinpi (double x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_d (y, x, MPFR_RNDN);
  int inex = mpfr_sinpi (y, y, rnd);
  mpfr_subnormalize (y, inex, rnd);
  double r = mpfr_get_d (y, MPFR_RNDN);
  mpfr_clear (y);
  return r;
}

double
ref_pow (double x, double y, mpfr_rnd_t rnd)
{
  mpfr_t z, _x, _y;
  int underflow = mpfr_flags_test (MPFR_FLAGS_UNDERFLOW);
  mpfr_inits2 (INTERNAL_PRECISION, z, _x, _y, NULL);
  mpfr_set_d (_x, x, MPFR_RNDN);
  mpfr_set_d (_y, y, MPFR_RNDN);
  int inex = mpfr_pow (z, _x, _y, rnd);
  inex = mpfr_subnormalize (z, inex, rnd);
  /* Workaround for bug in mpfr_subnormalize for MPFR <= 4.2.1: the underflow
     flag is set by mpfr_subnormalize() even for exact results. */
  if (inex == 0 && !underflow)
    mpfr_flags_clear (MPFR_FLAGS_UNDERFLOW);
  double ret = mpfr_get_d (z, rnd);
  mpfr_clears (z, _x, _y, NULL);
  return ret;
}

double
ref_tan (double x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
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
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_d (y, x, MPFR_RNDN);
  int inex = mpfr_tanh (y, y, rnd);
  mpfr_subnormalize (y, inex, rnd);
  double ret = mpfr_get_d (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

double
ref_tanpi (double x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_d (y, x, MPFR_RNDN);
  int inex = mpfr_tanpi (y, y, rnd);
  mpfr_subnormalize (y, inex, rnd);
  double ret = mpfr_get_d (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

double
ref_tgamma (double x, mpfr_rnd_t rnd)
{
  double fx = __builtin_floor (x);
  if (fx == x)
    if (x < 0.0)
      return __builtin_nanf ("12");
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_d (y, x, MPFR_RNDN);
  int inex = mpfr_gamma (y, y, rnd);
  mpfr_subnormalize (y, inex, rnd);
  double ret = mpfr_get_d (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}
