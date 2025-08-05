#include <float.h>
#include <math.h>
#include <stdint.h>
// NB: stdint should be included prior mpfr.h to define intmax_t, which
// mpfr.h uses to define mpfr_pown
#include <mpfr.h>

enum
{
  INTERNAL_PRECISION = FLT_MANT_DIG
};

float
ref_acosf (float x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_flt (y, x, MPFR_RNDN);
  int inex = mpfr_acos (y, y, rnd);
  mpfr_subnormalize (y, inex, rnd);
  float ret = mpfr_get_flt (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

float
ref_acoshf (float x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_flt (y, x, MPFR_RNDN);
  int inex = mpfr_acosh (y, y, rnd);
  mpfr_subnormalize (y, inex, rnd);
  float ret = mpfr_get_flt (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

float
ref_acospif (float x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_flt (y, x, MPFR_RNDN);
  mpfr_acospi (y, y, rnd);
  float ret = mpfr_get_flt (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

float
ref_asinf (float x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_flt (y, x, MPFR_RNDN);
  int inex = mpfr_asin (y, y, rnd);
  mpfr_subnormalize (y, inex, rnd);
  float ret = mpfr_get_flt (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

float
ref_asinhf (float x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_flt (y, x, MPFR_RNDN);
  int inex = mpfr_asinh (y, y, rnd);
  mpfr_subnormalize (y, inex, rnd);
  float ret = mpfr_get_flt (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

float
ref_asinpif (float x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_flt (y, x, MPFR_RNDN);
  int inex = mpfr_asinpi (y, y, rnd);
  mpfr_subnormalize (y, inex, rnd);
  float ret = mpfr_get_flt (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

float
ref_atanf (float x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_flt (y, x, MPFR_RNDN);
  int inex = mpfr_atan (y, y, rnd);
  mpfr_subnormalize (y, inex, rnd);
  float ret = mpfr_get_flt (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

float
ref_atan2f (float x, float y, mpfr_rnd_t rnd)
{
  mpfr_t xm, ym, zm;
  mpfr_init2 (xm, INTERNAL_PRECISION);
  mpfr_init2 (ym, INTERNAL_PRECISION);
  mpfr_init2 (zm, INTERNAL_PRECISION);
  mpfr_set_flt (xm, x, MPFR_RNDN);
  mpfr_set_flt (ym, y, MPFR_RNDN);
  int inex = mpfr_atan2 (zm, xm, ym, rnd);
  mpfr_subnormalize (zm, inex, rnd);
  float ret = mpfr_get_flt (zm, MPFR_RNDN);
  mpfr_clear (xm);
  mpfr_clear (ym);
  mpfr_clear (zm);
  return ret;
}

float
ref_atanhf (float x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_flt (y, x, MPFR_RNDN);
  int inex = mpfr_atanh (y, y, rnd);
  mpfr_subnormalize (y, inex, rnd);
  float ret = mpfr_get_flt (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

float
ref_atanpif (float x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_flt (y, x, MPFR_RNDN);
  int inex = mpfr_atanpi (y, y, rnd);
  mpfr_subnormalize (y, inex, rnd);
  float ret = mpfr_get_flt (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

float
ref_cbrtf (float x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_flt (y, x, MPFR_RNDN);
  int inex = mpfr_cbrt (y, y, rnd);
  mpfr_subnormalize (y, inex, rnd);
  float ret = mpfr_get_flt (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

float
ref_compoundnf (float x, long long int y, mpfr_rnd_t rnd)
{
  mpfr_t xm, zm;
  mpfr_init2 (xm, INTERNAL_PRECISION);
  mpfr_init2 (zm, INTERNAL_PRECISION);
  mpfr_set_flt (xm, x, MPFR_RNDN);
  int inex = mpfr_compound_si (zm, xm, y, rnd);
  mpfr_subnormalize (zm, inex, rnd);
  float ret = mpfr_get_flt (zm, MPFR_RNDN);
  mpfr_clear (xm);
  mpfr_clear (zm);
  return ret;
}

float
ref_cosf (float x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_flt (y, x, MPFR_RNDN);
  int inex = mpfr_cos (y, y, rnd);
  mpfr_subnormalize (y, inex, rnd);
  float ret = mpfr_get_flt (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

float
ref_coshf (float x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_flt (y, x, MPFR_RNDN);
  int inex = mpfr_cosh (y, y, rnd);
  mpfr_subnormalize (y, inex, rnd);
  float ret = mpfr_get_flt (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

float
ref_cospif (float x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_flt (y, x, MPFR_RNDN);
  mpfr_cospi (y, y, rnd);
  float ret = mpfr_get_flt (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

float
ref_erff (float x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_flt (y, x, MPFR_RNDN);
  int inex = mpfr_erf (y, y, rnd);
  mpfr_subnormalize (y, inex, rnd);
  float ret = mpfr_get_flt (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

float
ref_erfcf (float x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_flt (y, x, MPFR_RNDN);
  int inex = mpfr_erfc (y, y, rnd);
  mpfr_subnormalize (y, inex, rnd);
  float ret = mpfr_get_flt (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

float
ref_exp10f (float x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_flt (y, x, MPFR_RNDN);
  int inex = mpfr_exp10 (y, y, rnd);
  mpfr_subnormalize (y, inex, rnd);
  float ret = mpfr_get_flt (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

float
ref_exp10m1f (float x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_flt (y, x, MPFR_RNDN);
  int inex = mpfr_exp10m1 (y, y, rnd);
  mpfr_subnormalize (y, inex, rnd);
  float ret = mpfr_get_flt (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

float
ref_exp2f (float x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_flt (y, x, MPFR_RNDN);
  int inex = mpfr_exp2 (y, y, rnd);
  mpfr_subnormalize (y, inex, rnd);
  float ret = mpfr_get_flt (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

float
ref_exp2m1f (float x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_flt (y, x, MPFR_RNDN);
  int inex = mpfr_exp2m1 (y, y, rnd);
  mpfr_subnormalize (y, inex, rnd);
  float ret = mpfr_get_flt (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

float
ref_expf (float x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_flt (y, x, MPFR_RNDN);
  int inex = mpfr_exp (y, y, rnd);
  mpfr_subnormalize (y, inex, rnd);
  float ret = mpfr_get_flt (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

float
ref_expm1f (float x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_flt (y, x, MPFR_RNDN);
  int inex = mpfr_expm1 (y, y, rnd);
  mpfr_subnormalize (y, inex, rnd);
  float ret = mpfr_get_flt (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

float
ref_hypotf (float x, float y, mpfr_rnd_t rnd)
{
  mpfr_t xm, ym, zm;

  union
  {
    float f;
    uint32_t u;
  } xi = { .f = x }, yi = { .f = y };
  if ((xi.u << 1) < (0xff8ull << 20)
      && (xi.u << 1) > (0xff0ull << 20)) // x = sNAN
    return x + y;                        // will return qNaN
  if ((yi.u << 1) < (0xff8ull << 20)
      && (yi.u << 1) > (0xff0ull << 20)) // y = sNAN
    return x + y;                        // will return qNaN
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

  mpfr_init2 (xm, INTERNAL_PRECISION);
  mpfr_init2 (ym, INTERNAL_PRECISION);
  mpfr_init2 (zm, INTERNAL_PRECISION);
  mpfr_set_flt (xm, x, MPFR_RNDN);
  mpfr_set_flt (ym, y, MPFR_RNDN);
  int inex = mpfr_hypot (zm, xm, ym, rnd);
  mpfr_subnormalize (zm, inex, rnd);
  float ret = mpfr_get_flt (zm, MPFR_RNDN);
  mpfr_clear (xm);
  mpfr_clear (ym);
  mpfr_clear (zm);
  return ret;
}

float
ref_lgammaf (float x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_flt (y, x, MPFR_RNDN);
  int inex = mpfr_lgamma (y, &signgam, y, rnd);
  mpfr_subnormalize (y, inex, rnd);
  float ret = mpfr_get_flt (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

float
ref_logf (float x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_flt (y, x, MPFR_RNDN);
  int inex = mpfr_log (y, y, rnd);
  mpfr_subnormalize (y, inex, rnd);
  float ret = mpfr_get_flt (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

float
ref_log1pf (float x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_flt (y, x, MPFR_RNDN);
  int inex = mpfr_log1p (y, y, rnd);
  mpfr_subnormalize (y, inex, rnd);
  float ret = mpfr_get_flt (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

float
ref_log2f (float x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_flt (y, x, MPFR_RNDN);
  int inex = mpfr_log2 (y, y, rnd);
  mpfr_subnormalize (y, inex, rnd);
  float ret = mpfr_get_flt (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

float
ref_log2p1f (float x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_flt (y, x, MPFR_RNDN);
  int inex = mpfr_log2p1 (y, y, rnd);
  mpfr_subnormalize (y, inex, rnd);
  float ret = mpfr_get_flt (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

float
ref_log10f (float x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_flt (y, x, MPFR_RNDN);
  int inex = mpfr_log10 (y, y, rnd);
  mpfr_subnormalize (y, inex, rnd);
  float ret = mpfr_get_flt (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

float
ref_log10p1f (float x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_flt (y, x, MPFR_RNDN);
  int inex = mpfr_log10p1 (y, y, rnd);
  mpfr_subnormalize (y, inex, rnd);
  float ret = mpfr_get_flt (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

float
ref_rsqrtf (float x, mpfr_rnd_t rnd)
{
  /* mpfr_rec_sqrt differs from IEEE 754-2019: IEEE 754-2019 says that
     rsqrt(-0) should give -Inf, whereas mpfr_rec_sqrt(-0) gives +Inf */
  if (x == 0.0f && 1.0f / x < 0.0f)
    return 1.0f / x;
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_flt (y, x, MPFR_RNDN);
  mpfr_rec_sqrt (y, y, rnd);
  /* since |x| < 2^128 for non-zero x, we have 2^-64 < 1/sqrt(x),
     thus no underflow can happen, and there is no need to call
     mpfr_subnormalize */
  float ret = mpfr_get_flt (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

float
ref_sinf (float x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_flt (y, x, MPFR_RNDN);
  int inex = mpfr_sin (y, y, rnd);
  mpfr_subnormalize (y, inex, rnd);
  float ret = mpfr_get_flt (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

float
ref_sinhf (float x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_flt (y, x, MPFR_RNDN);
  int inex = mpfr_sinh (y, y, rnd);
  mpfr_subnormalize (y, inex, rnd);
  float ret = mpfr_get_flt (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

float
ref_sinpif (float x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_flt (y, x, MPFR_RNDN);
  int inex = mpfr_sinpi (y, y, rnd);
  mpfr_subnormalize (y, inex, rnd);
  float ret = mpfr_get_flt (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

float
ref_powf (float x, float y, mpfr_rnd_t rnd)
{
  mpfr_t xm, ym, zm;
  mpfr_init2 (xm, INTERNAL_PRECISION);
  mpfr_init2 (ym, INTERNAL_PRECISION);
  mpfr_init2 (zm, INTERNAL_PRECISION);
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

float
ref_pownf (float x, long long int y, mpfr_rnd_t rnd)
{
  mpfr_t xm, zm;
  mpfr_init2 (xm, INTERNAL_PRECISION);
  mpfr_init2 (zm, INTERNAL_PRECISION);
  mpfr_set_flt (xm, x, MPFR_RNDN);
  int inex = mpfr_pown (zm, xm, y, rnd);
  mpfr_subnormalize (zm, inex, rnd);
  float ret = mpfr_get_flt (zm, MPFR_RNDN);
  mpfr_clear (xm);
  mpfr_clear (zm);
  return ret;
}

float
ref_tanf (float x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_flt (y, x, MPFR_RNDN);
  int inex = mpfr_tan (y, y, rnd);
  mpfr_subnormalize (y, inex, rnd);
  float ret = mpfr_get_flt (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

float
ref_tanhf (float x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_flt (y, x, MPFR_RNDN);
  int inex = mpfr_tanh (y, y, rnd);
  mpfr_subnormalize (y, inex, rnd);
  float ret = mpfr_get_flt (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

float
ref_tanpif (float x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_flt (y, x, MPFR_RNDN);
  int inex = mpfr_tanpi (y, y, rnd);
  mpfr_subnormalize (y, inex, rnd);
  float ret = mpfr_get_flt (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

float
ref_tgammaf (float x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_init2 (y, INTERNAL_PRECISION);
  mpfr_set_flt (y, x, MPFR_RNDN);
  int inex = mpfr_gamma (y, y, rnd);
  mpfr_subnormalize (y, inex, rnd);
  float ret = mpfr_get_flt (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}
