#ifndef REFIMPLS_H
#define REFIMPLS_H

#include <expected>
#include <fenv.h>
#include <mpfr.h>
#include <string_view>

#include "cxxcompat.h"

namespace refimpls
{

template <typename F> using func_f_t = F (*) (F);
template <typename F> using func_f_mpfr_t = F (*) (F, mpfr_rnd_t);

template <typename F> using func_f_f_t = F (*) (F, F);
template <typename F> using func_f_f_mpfr_t = F (*) (F, F, mpfr_rnd_t);

template <typename F> using func_f_lli_t = F (*) (F, long long int);
template <typename F>
using func_f_lli_mpfr_t = F (*) (F, long long int, mpfr_rnd_t);

template <typename T> struct func_ref_t
{
  func_ref_t (const func_f_mpfr_t<T> &func) : f (func) {}

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

  const func_f_mpfr_t<T> f;
};

template <typename T> struct func_f_f_ref_t
{
  func_f_f_ref_t (const func_f_f_mpfr_t<T> &func) : f (func) {}

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

  const func_f_f_mpfr_t<T> f;
};

template <typename T> struct func_f_lli_ref_t
{
  func_f_lli_ref_t (const func_f_lli_mpfr_t<T> &func) : f (func) {}

  T
  operator() (T x, long long int y, int rnd) const
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

  const func_f_lli_mpfr_t<T> f;
};

enum class errors_t
{
  invalid_func
};

enum class func_type_t
{
  f32_f,
  f32_f_f,
  f32_f_lli,
  f64_f,
  f64_f_f,
  f64_f_lli,
};

template <class F> void setup_ref_impl ();

std::expected<func_type_t, errors_t> get_func_type (const std::string_view &);

template <typename F>
std::expected<std::pair<func_f_t<F>, func_ref_t<F> >, errors_t>
get_f (const std::string_view &);

template <typename F>
std::expected<std::pair<func_f_f_t<F>, func_f_f_ref_t<F> >, errors_t>
get_f_f (const std::string_view &);

template <typename F>
std::expected<std::pair<func_f_lli_t<F>, func_f_lli_ref_t<F> >, errors_t>
get_f_lli (const std::string_view &);

} // refimpls

#endif
