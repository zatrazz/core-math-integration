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
#include "refimpls_modes.h"

namespace refimpls
{

// The reference (MPFR) functions compute the multi-precision result once and
// round it to every rounding mode selected in a bitmask (1u << REF_RND*),
// storing each rounded result in out[REF_RND*].  See refimpls_modes.h.

template <typename F> using FuncF = F (*) (F);
template <typename F> using FuncFMpfr = void (*) (F, unsigned, F[REF_NRND]);

template <typename F> using FuncFpFp = void (*) (F, F *, F *);
template <typename F>
using FuncFpFpMpfr = void (*) (F, unsigned, F[REF_NRND], F[REF_NRND]);

template <typename F> using FuncFF = F (*) (F, F);
template <typename F> using FuncFFMpfr = void (*) (F, F, unsigned, F[REF_NRND]);

template <typename F> using FuncFLLI = F (*) (F, long long int);
template <typename F>
using FuncFLLIMpfr = void (*) (F, long long int, unsigned, F[REF_NRND]);

template <typename T> struct FuncFReference
{
  FuncFReference (const FuncFMpfr<T> &func) : f (func) {}

  void
  operator() (T input, unsigned mask, T out[REF_NRND]) const
  {
    f (input, mask, out);
  }

  const FuncFMpfr<T> f;
};

template <typename T> struct FuncFFReference
{
  FuncFFReference (const FuncFFMpfr<T> &func) : f (func) {}

  void
  operator() (T x, T y, unsigned mask, T out[REF_NRND]) const
  {
    f (x, y, mask, out);
  }

  const FuncFFMpfr<T> f;
};

template <typename T> struct FuncFpFpReference
{
  FuncFpFpReference (const FuncFpFpMpfr<T> &func) : f (func) {}

  void
  operator() (T x, unsigned mask, T out1[REF_NRND], T out2[REF_NRND]) const
  {
    f (x, mask, out1, out2);
  }

  const FuncFpFpMpfr<T> f;
};

template <typename T> struct FuncFLLIReference
{
  FuncFLLIReference (const FuncFLLIMpfr<T> &func) : f (func) {}

  void
  operator() (T x, long long int y, unsigned mask, T out[REF_NRND]) const
  {
    f (x, y, mask, out);
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
