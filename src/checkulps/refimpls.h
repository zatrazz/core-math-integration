#ifndef REFIMPLS_H
#define REFIMPLS_H

#include <expected>
#include <mpfr.h>
#include <string_view>

#include "cxxcompat.h"

namespace refimpls
{

template <typename F> using univariate_t = F (*) (F);
template <typename F> using univariate_mpfr_t = F (*) (F, mpfr_rnd_t);

template <typename F> using bivariate_t = F (*) (F, F);
template <typename F> using bivariate_mpfr_t = F (*) (F, F, mpfr_rnd_t);

template <typename T> struct univariate_ref_t
{
  univariate_ref_t (const univariate_mpfr_t<T> &func) : f (func) {}

  double
  operator() (T input, int rnd) const
  {
    switch (rnd)
      {
      case FE_TONEAREST:
        return f (input, MPFR_RNDN);
      case FE_UPWARD:
        return f (input, MPFR_RNDU);
      case FE_DOWNWARD:
        return f (input, MPFR_RNDD);
      case FE_TOWARDZERO:
        return f (input, MPFR_RNDZ);
      default:
        std::unreachable ();
      };
  }

  const univariate_mpfr_t<T> f;
};

template <typename T> struct bivariate_ref_t
{
  bivariate_ref_t (const bivariate_mpfr_t<T> &func) : f (func) {}

  T
  operator() (T x, T y, int rnd) const
  {
    switch (rnd)
      {
      case FE_TONEAREST:
        return f (x, y, MPFR_RNDN);
      case FE_UPWARD:
        return f (x, y, MPFR_RNDU);
      case FE_DOWNWARD:
        return f (x, y, MPFR_RNDD);
      case FE_TOWARDZERO:
        return f (x, y, MPFR_RNDZ);
      default:
        std::unreachable ();
      };
  }

  const bivariate_mpfr_t<T> f;
};

enum class errors_t
{
  invalid_func
};

enum class func_type_t
{
  binary32_univariate,
  binary32_bivariate,
  binary64_univariate,
  binary64_bivariate,
};

template <class F> void setup_ref_impl ();

std::expected<func_type_t, errors_t> get_func_type (const std::string_view &);

template <typename F>
std::expected<std::pair<univariate_t<F>, univariate_ref_t<F> >, errors_t>
get_univariate (const std::string_view &, bool coremath = false);

template <typename F>
std::expected<std::pair<bivariate_t<F>, bivariate_ref_t<F> >, errors_t>
get_bivariate (const std::string_view &, bool coremath = false);

} // refimpls

#endif
