#include <utility>
#include <cmath>
#include <cstdint>

#include <fenv.h>
#include <mpfr.h>

#include "refimpls.h"

namespace refimpls
{

extern "C" {
  // From core-math
  double cr_asin (double);

  double acospi (double) __attribute__ ((weak)); 
  double asinpi (double) __attribute__ ((weak)); 
};;

template <auto F, typename T>
class ref_mode_univariate
{
public:
  static T tonearest  (T x) { return F(x, MPFR_RNDN); }
  static T toupward   (T x) { return F(x, MPFR_RNDU); }
  static T todownward (T x) { return F(x, MPFR_RNDD); }
  static T towardzero (T x) { return F(x, MPFR_RNDZ); }

  static univariate_t get (int rnd)
  {
    switch (rnd)
    {
    case FE_TONEAREST:  return ref_mode_univariate<F, T>::tonearest;
    case FE_UPWARD:     return ref_mode_univariate<F, T>::toupward;
    case FE_DOWNWARD:   return ref_mode_univariate<F, T>::todownward;
    case FE_TOWARDZERO: return ref_mode_univariate<F, T>::towardzero;
    default:	        std::unreachable();
    }
  }
};

static double
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

typedef union {double f; std::uint64_t u;} b64u64_u;

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
#if 0
  if (std::isnan (x))
    return x;
#endif
  mpfr_t y;
  mpfr_init2(y, 53);
  mpfr_set_d(y, x, MPFR_RNDN);
  int d = mpfr_asinh(y, y, rnd);
  mpfr_subnormalize(y, d, rnd);
  double ret = mpfr_get_d(y, MPFR_RNDN);
  mpfr_clear(y);
  return ret;
}

std::expected<univariate_t, errors_t>
get_univariate (const std::string_view &str, bool coremath)
{
  if (str == "acos")
    return acos;
  else if (str == "acosh")
    return acosh;
  else if (str == "asin")
    return coremath ? cr_asin : asin;
  else if (str == "asinh")
    return asinh;
  else if (str == "acospi" && acospi != nullptr)
    return acospi;
  else if (str == "asinpi" && asinpi != nullptr)
    return asinpi;
  return std::unexpected (errors_t::invalid_func);
}

std::expected<univariate_t, errors_t>
get_univariate_ref (const std::string_view &str, int rnd)
{
  if (str == "acos")
    return ref_mode_univariate<ref_acos, double>::get(rnd);
  if (str == "acospi")
    return ref_mode_univariate<ref_acospi, double>::get(rnd);
  else if (str == "acosh")
    return ref_mode_univariate<ref_acosh, double>::get(rnd);
  else if (str == "asin")
    return ref_mode_univariate<ref_asin, double>::get(rnd);
  else if (str == "asinpi")
    return ref_mode_univariate<ref_asinpi, double>::get(rnd);
  else if (str == "asinh")
    return ref_mode_univariate<ref_asinh, double>::get(rnd);
  return std::unexpected (errors_t::invalid_func);
}

}
