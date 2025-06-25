#include <utility>
#include <cmath>

#include <fenv.h>
#include <mpfr.h>

#include "refimpls.h"

namespace refimpls
{

static inline mpfr_rnd_t
randmode_to_mpfr (int rand)
{
  switch (rand)
  {
  case FE_TONEAREST:  return MPFR_RNDN;
  case FE_UPWARD:     return MPFR_RNDU;
  case FE_DOWNWARD:   return MPFR_RNDD;
  case FE_TOWARDZERO: return MPFR_RNDZ;
  default:	      std::unreachable();
  }
}

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
ref_acos (double x, int rnd)
{
  mpfr_rnd_t r = randmode_to_mpfr (rnd);
  mpfr_t y;
  mpfr_init2 (y, 53);
  mpfr_set_d (y, x, MPFR_RNDN);
  int inex = mpfr_acos (y, y, r);
  mpfr_subnormalize (y, inex, r);
  double ret = mpfr_get_d (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

std::expected<univariate_t, errors_t>
get_univariate (const std::string_view &str)
{
  if (str == "acos")
    return acos;
  return std::unexpected (errors_t::invalid_func);
}

std::expected<univariate_t, errors_t>
get_univariate_ref (const std::string_view &str, int rnd)
{
  if (str == "acos")
    return ref_mode_univariate<ref_acos, double>::get(rnd);
  return std::unexpected (errors_t::invalid_func);
}

}
