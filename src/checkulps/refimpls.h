//
// Copyright (c) Adhemerval Zanella. All rights reserved.
// Licensed under the MIT license. See LICENSE file in the project root for
// details.
//

#ifndef REFIMPLS_H
#define REFIMPLS_H

#include <expected>
#include <fenv.h>
#include <mpfr.h>
#include <format>
#include <string_view>

#include "cxxcompat.h"

namespace refimpls
{

template <typename F> using func_f_t = F (*) (F);
template <typename F> using func_f_mpfr_t = F (*) (F, mpfr_rnd_t);

template <typename F> using func_f_fp_fp_t = void (*) (F, F *, F *);
template <typename F>
using func_f_fp_fp_mpfr_t = void (*) (F, F *, F *, mpfr_rnd_t);

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

template <typename T> struct func_f_fp_fp_ref_t
{
  func_f_fp_fp_ref_t (const func_f_fp_fp_mpfr_t<T> &func) : f (func) {}

  void
  operator() (T x, T *r1, T *r2, int rnd) const
  {
    switch (rnd)
      {
      case FE_TONEAREST:
	return f (x, r1, r2, MPFR_RNDN);
      case FE_UPWARD:
	return f (x, r1, r2, MPFR_RNDU);
      case FE_DOWNWARD:
	return f (x, r1, r2, MPFR_RNDD);
      case FE_TOWARDZERO:
	return f (x, r1, r2, MPFR_RNDZ);
      default:
	std::unreachable ();
      };
  }

  const func_f_fp_fp_mpfr_t<T> f;
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
  f32_f_fp_fp,
  f32_f_lli,
  f64_f,
  f64_f_f,
  f64_f_fp_fp,
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
std::expected<std::pair<func_f_fp_fp_t<F>, func_f_fp_fp_ref_t<F> >, errors_t>
get_f_fp_fp (const std::string_view &);

template <typename F>
std::expected<std::pair<func_f_lli_t<F>, func_f_lli_ref_t<F> >, errors_t>
get_f_lli (const std::string_view &);

} // refimpls

template <>
struct std::formatter<refimpls::func_type_t> : std::formatter<std::string_view>
{
  auto
  format (refimpls::func_type_t t, std::format_context &ctx) const
  {
    std::string_view r;
    switch (t)
      {
      case refimpls::func_type_t::f32_f:
	r = "float (*)(float)";
	break;
      case refimpls::func_type_t::f32_f_f:
	r = "float (*)(float, float)";
	break;
      case refimpls::func_type_t::f32_f_fp_fp:
	r = "void (*)(float, float*, float*)";
	break;
      case refimpls::func_type_t::f32_f_lli:
	r = "float (*)(float, long long)";
	break;
      case refimpls::func_type_t::f64_f:
	r = "double (*)(double)";
	break;
      case refimpls::func_type_t::f64_f_f:
	r = "double (*)(double, double)";
	break;
      case refimpls::func_type_t::f64_f_fp_fp:
	r = "void (*)(double, double*, double*)";
	break;
      case refimpls::func_type_t::f64_f_lli:
	r = "double (*)(double, long long)";
	break;
      }
    return std::formatter<std::string_view>::format (r, ctx);
  }
};

#endif
