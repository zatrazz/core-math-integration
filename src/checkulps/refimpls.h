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

template <typename F> using FuncF = F (*) (F);
template <typename F> using FuncFMpfr = F (*) (F, mpfr_rnd_t);

template <typename F> using FuncFpFp = void (*) (F, F *, F *);
template <typename F> using FuncFpFpMpfr = void (*) (F, F *, F *, mpfr_rnd_t);

template <typename F> using FuncFF = F (*) (F, F);
template <typename F> using FuncFFMpfr = F (*) (F, F, mpfr_rnd_t);

template <typename F> using FuncFLLI = F (*) (F, long long int);
template <typename F>
using FuncFLLIMpfr = F (*) (F, long long int, mpfr_rnd_t);

template <typename T> struct FuncFReference
{
  FuncFReference (const FuncFMpfr<T> &func) : f (func) {}

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

  const FuncFMpfr<T> f;
};

template <typename T> struct FuncFFReference
{
  FuncFFReference (const FuncFFMpfr<T> &func) : f (func) {}

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

  const FuncFFMpfr<T> f;
};

template <typename T> struct FuncFpFpReference
{
  FuncFpFpReference (const FuncFpFpMpfr<T> &func) : f (func) {}

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

  const FuncFpFpMpfr<T> f;
};

template <typename T> struct FuncFLLIReference
{
  FuncFLLIReference (const FuncFLLIMpfr<T> &func) : f (func) {}

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

  const FuncFLLIMpfr<T> f;
};

enum class Errors
{
  invalid_func
};

enum class FunctionType
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

template <class F> void setupReferenceImpl ();

std::expected<FunctionType, Errors> getFunctionType (const std::string_view &);

template <typename F>
std::expected<std::pair<FuncF<F>, FuncFReference<F> >, Errors>
getFunctionFloat (const std::string_view &);

template <typename F>
std::expected<std::pair<FuncFF<F>, FuncFFReference<F> >, Errors>
getFunctionFloatFloat (const std::string_view &);

template <typename F>
std::expected<std::pair<FuncFpFp<F>, FuncFpFpReference<F> >, Errors>
getFunctionFloatpFloatp (const std::string_view &);

template <typename F>
std::expected<std::pair<FuncFLLI<F>, FuncFLLIReference<F> >, Errors>
getFunctionFloatLLI (const std::string_view &);

} // refimpls

template <>
struct std::formatter<refimpls::FunctionType>
    : std::formatter<std::string_view>
{
  auto
  format (refimpls::FunctionType t, std::format_context &ctx) const
  {
    std::string_view r;
    switch (t)
      {
      case refimpls::FunctionType::f32_f:
	r = "float (*)(float)";
	break;
      case refimpls::FunctionType::f32_f_f:
	r = "float (*)(float, float)";
	break;
      case refimpls::FunctionType::f32_f_fp_fp:
	r = "void (*)(float, float*, float*)";
	break;
      case refimpls::FunctionType::f32_f_lli:
	r = "float (*)(float, long long)";
	break;
      case refimpls::FunctionType::f64_f:
	r = "double (*)(double)";
	break;
      case refimpls::FunctionType::f64_f_f:
	r = "double (*)(double, double)";
	break;
      case refimpls::FunctionType::f64_f_fp_fp:
	r = "void (*)(double, double*, double*)";
	break;
      case refimpls::FunctionType::f64_f_lli:
	r = "double (*)(double, long long)";
	break;
      }
    return std::formatter<std::string_view>::format (r, ctx);
  }
};

#endif
