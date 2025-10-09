//
// Copyright (c) Adhemerval Zanella. All rights reserved.
// Licensed under the MIT license. See LICENSE file in the project root for
// details.
//

#include <algorithm>
#include <array>
#include <cmath>
#include <string>
#include <utility>

#include <fenv.h>

#include "config.h"
#include "refimpls.h"

#ifndef HAVE_LGAMMAF_R_PROTOTYE
extern "C"
{
  float lgammaf_r (float x, int *psigngam);
};
#endif
#ifndef HAVE_LGAMMA_R_PROTOTYE
extern "C"
{
  double lgamma_r (double x, int *psigngam);
};
#endif

namespace refimpls
{

template <>
void
setupReferenceImpl<float> ()
{
  mpfr_set_emin (-148);
  mpfr_set_emax (128);
}

template <>
void
setupReferenceImpl<double> ()
{
  mpfr_set_emin (-1073);
  mpfr_set_emax (1024);
}

extern "C"
{
#define _DEF_F(__name)                                                        \
  extern float ref_##__name##f (float, mpfr_rnd_t);                           \
  extern double ref_##__name (double, mpfr_rnd_t)

#define DEF_F(__name) _DEF_F (__name)

#define _DEF_F_FP_FP(__name)                                                  \
  extern void ref_##__name##f (float, float *, float *, mpfr_rnd_t)           \
      __attribute__ ((weak));                                                 \
  extern void ref_##__name (double, double *, double *, mpfr_rnd_t)           \
      __attribute__ ((weak))

#define DEF_F_FP_FP_WEAK(__name)                                              \
  extern void __name##f (float, float *, float *)                             \
      __attribute__ ((weak, used));                                           \
  extern void __name (double, double *, double *)                             \
      __attribute__ ((weak, used));                                           \
  _DEF_F_FP_FP (__name)

#define DEF_F_WEAK(__name)                                                    \
  extern float __name##f (float) __attribute__ ((weak));                      \
  extern double __name (double) __attribute__ ((weak));                       \
  _DEF_F (__name)

#define _DEF_F_F(__name)                                                      \
  extern float ref_##__name##f (float, float, mpfr_rnd_t);                    \
  extern double ref_##__name (double, double, mpfr_rnd_t)

#define DEF_F_F(__name) _DEF_F_F (__name)

#define DEF_F_F_WEAK(__name)                                                  \
  extern float __name##f (float, float) __attribute__ ((weak, used));         \
  extern double __name (double, double) __attribute__ ((weak, used));         \
  _DEF_F_F (__name)

#define _DEF_F_LI(__name)                                                     \
  extern float ref_##__name##f (float, long long int, mpfr_rnd_t);            \
  extern double ref_##__name (double, long long int, mpfr_rnd_t)

#define DEF_F_LLI_WEAK(__name)                                                \
  extern float __name##f (float, long long int) __attribute__ ((weak, used)); \
  extern double __name (double, long long int) __attribute__ ((weak, used));  \
  _DEF_F_LI (__name)

  DEF_F_F (atan2);
  DEF_F_WEAK (atanpi);
  DEF_F (acos);
  DEF_F (acosh);
  DEF_F_WEAK (acospi);
  DEF_F (asin);
  DEF_F (asinh);
  DEF_F_WEAK (asinpi);
  DEF_F (atan);
  DEF_F (atanh);
  DEF_F (cbrt);
  DEF_F_LLI_WEAK (compoundn);
  DEF_F (cos);
  DEF_F (cosh);
  DEF_F_WEAK (cospi);
  DEF_F (erf);
  DEF_F (erfc);
  DEF_F (exp);
  DEF_F (expm1);
  DEF_F_WEAK (exp10);
  DEF_F_WEAK (exp10m1);
  DEF_F (exp2);
  DEF_F_WEAK (exp2m1);
  DEF_F (lgamma);
  DEF_F (log);
  DEF_F (log1p);
  DEF_F (log2);
  DEF_F_WEAK (log2p1);
  DEF_F (log10);
  DEF_F_WEAK (log10p1);
  DEF_F_F (hypot);
  DEF_F_F (pow);
  DEF_F_LLI_WEAK (pown);
  DEF_F_F_WEAK (powr);
  DEF_F_LLI_WEAK (rootn);
  DEF_F_WEAK (rsqrt);
  DEF_F (sin);
  DEF_F_FP_FP_WEAK (sincos);
  DEF_F (sinh);
  DEF_F_WEAK (sinpi);
  DEF_F (tan);
  DEF_F (tanh);
  DEF_F_WEAK (tanpi);
  DEF_F (tgamma);

#undef _DEF_F
#undef DEF_F
#undef DEF_F_WEAK
#undef _DEF_F_F
#undef DEF_F_F
#undef DEF_F_F_WEAK
};

template <typename F, typename F_MPFR> struct FuncFDescription
{
  const std::string name;
  volatile F func;
  F_MPFR mpfrFunc;
};

typedef FuncFDescription<FuncF<float>, FuncFMpfr<float> > FuncF32Description;

static float
lgammaf_wrapper (float x)
{
  // The lgamma is not thread-safe, and we can ignore the sign.
  int sign;
  return lgammaf_r (x, &sign);
}

// clang-format off
const static std::array funcF32 = {
#define FUNC_DEF(name)        \
  FuncF32Description {    \
    #name, name, ref_##name   \
  }
  FUNC_DEF (atanpif),
  FUNC_DEF (acosf),
  FUNC_DEF (acoshf),
  FUNC_DEF (acospif),
  FUNC_DEF (asinf),
  FUNC_DEF (asinhf),
  FUNC_DEF (asinpif),
  FUNC_DEF (atanf),
  FUNC_DEF (atanhf),
  FUNC_DEF (cbrtf),
  FUNC_DEF (cosf),
  FUNC_DEF (coshf),
  FUNC_DEF (cospif),
  FUNC_DEF (erff),
  FUNC_DEF (erfcf),
  FUNC_DEF (expf),
  FUNC_DEF (expm1f),
  FUNC_DEF (exp10f),
  FUNC_DEF (exp10m1f),
  FUNC_DEF (exp2f),
  FUNC_DEF (exp2m1f),
  FuncF32Description {
    "lgammaf", lgammaf_wrapper, ref_lgammaf
  },
  FUNC_DEF (logf),
  FUNC_DEF (log1pf),
  FUNC_DEF (log2f),
  FUNC_DEF (log2p1f),
  FUNC_DEF (log10f),
  FUNC_DEF (log10p1f),
  FUNC_DEF (rsqrtf),
  FUNC_DEF (sinf),
  FUNC_DEF (sinhf),
  FUNC_DEF (sinpif),
  FUNC_DEF (tanf),
  FUNC_DEF (tanhf),
  FUNC_DEF (tanpif),
  FUNC_DEF (tgammaf),
#undef FUNC_DEF
  };
// clang-format on

typedef FuncFDescription<FuncF<double>, FuncFMpfr<double> > FuncF64Description;

static double
lgamma_wrapper (double x)
{
  // The lgamma is not thread-safe, and we can ignore the sign.
  int sign;
  return lgamma_r (x, &sign);
}

// clang-format off
const static std::array funcF64 = {
#define FUNC_DEF(name)        \
  FuncF64Description {    \
    #name, name, ref_##name   \
  }
  FUNC_DEF (atanpi),
  FUNC_DEF (acos),
  FUNC_DEF (acosh),
  FUNC_DEF (acospi),
  FUNC_DEF (asin),
  FUNC_DEF (asinh),
  FUNC_DEF (asinpi),
  FUNC_DEF (atan),
  FUNC_DEF (atanh),
  FUNC_DEF (cbrt),
  FUNC_DEF (cos),
  FUNC_DEF (cosh),
  FUNC_DEF (cospi),
  FUNC_DEF (erf),
  FUNC_DEF (erfc),
  FUNC_DEF (exp),
  FUNC_DEF (expm1),
  FUNC_DEF (exp10),
  FUNC_DEF (exp10m1),
  FUNC_DEF (exp2),
  FUNC_DEF (exp2m1),
  FuncF64Description {
    "lgamma", lgamma_wrapper, ref_lgamma
  },
  FUNC_DEF (log),
  FUNC_DEF (log1p),
  FUNC_DEF (log2),
  FUNC_DEF (log2p1),
  FUNC_DEF (log10),
  FUNC_DEF (log10p1),
  FUNC_DEF (rsqrt),
  FUNC_DEF (sin),
  FUNC_DEF (sinh),
  FUNC_DEF (sinpi),
  FUNC_DEF (tan),
  FUNC_DEF (tanh),
  FUNC_DEF (tanpi),
  FUNC_DEF (tgamma),
#undef FUNC_DEF
  };
// clang-format on

typedef FuncFDescription<FuncFpFp<float>, FuncFpFpMpfr<float> >
    FuncFp32Fp32Description;

// clang-format off
const static std::array funcFP32FP32 = {
#define FUNC_DEF(name)                                               \
  FuncFp32Fp32Description {                                          \
    #name, name, ref_##name                                          \
  }
  FUNC_DEF (sincosf),
#undef FUNC_DEF
};
// clang-format on

typedef FuncFDescription<FuncFpFp<double>, FuncFpFpMpfr<double> >
    FuncFp64Fp64Description;

// clang-format off
const static std::array funcFP64FP64 = {
#define FUNC_DEF(name)                                               \
  FuncFp64Fp64Description {                                          \
    #name, name, ref_##name                                          \
  }
  FUNC_DEF (sincos),
#undef FUNC_DEF
};
// clang-format on

typedef FuncFDescription<FuncFF<float>, FuncFFMpfr<float> >
    FuncF32F32Description;

// clang-format off
const static std::array funcF32F32 = {
#define FUNC_DEF(name)                                              \
  FuncF32F32Description {                                           \
    #name, name, ref_##name                                         \
  }
  FUNC_DEF (atan2f),
  FUNC_DEF (hypotf),
  FUNC_DEF (powf),
  FUNC_DEF (powrf),
#undef FUNC_DEF
  };
// clang-format on

typedef FuncFDescription<FuncFF<double>, FuncFFMpfr<double> >
    FuncF64F64Description;

// clang-format off
const static std::array funcF64F64 = {
#define FUNC_DEF(name)                                             \
  FuncF64F64Description {                                          \
    #name, name, ref_##name                                        \
  }
  FUNC_DEF (atan2),
  FUNC_DEF (hypot),
  FUNC_DEF (pow),
  FUNC_DEF (powr),
#undef FUNC_DEF
  };
// clang-format on

typedef FuncFDescription<FuncFLLI<float>, FuncFLLIMpfr<float> >
    FuncF32LLIDescription;

// clang-format off
const static std::array funcF32LLI = {
#define FUNC_DEF(name)                                               \
  FuncF32LLIDescription {                                            \
    #name, name, ref_##name                                          \
  }
  FUNC_DEF (compoundnf),
  FUNC_DEF (pownf),
  FUNC_DEF (rootnf),
#undef FUNC_DEF
  };
// clang-format on

typedef FuncFDescription<FuncFLLI<double>, FuncFLLIMpfr<double> >
    FuncF64LLIDescription;

// clang-format off
const static std::array funcF64LLI = {
#define FUNC_DEF(name)                                               \
  FuncF64LLIDescription {                                            \
    #name, name, ref_##name                                          \
  }
  FUNC_DEF (compoundn),
  FUNC_DEF (pown),
  FUNC_DEF (rootn),
#undef FUNC_DEF
  };
// clang-format on

enum class FuncError
{
  FUNCTION_NOT_FOUND
};

template <typename T, std::size_t N>
static std::expected<typename std::array<T, N>::const_iterator, FuncError>
find_function (const std::array<T, N> &funcs, const std::string_view &funcname)
{
  auto it = std::find_if (
      funcs.begin (), funcs.end (),
      [&funcname] (const std::array<T, N>::const_reference &func) {
	return func.name == funcname;
      });
  if (it != funcs.end ())
    return it;
  return std::unexpected (FuncError::FUNCTION_NOT_FOUND);
}

template <typename T, std::size_t N>
static bool
contains_function (const std::array<T, N> &funcs,
		   const std::string_view &funcname)
{
  if (auto it = find_function (funcs, funcname))
    return true;
  return false;
}

template <>
std::expected<std::pair<FuncF<float>, FuncFReference<float> >, Errors>
getFunctionFloat (const std::string_view &funcname)
{
  if (const auto it = find_function (funcF32, funcname))
    return std::make_pair (*it.value ()->func,
			   FuncFReference<float>{ *it.value ()->mpfrFunc });
  return std::unexpected (Errors::invalid_func);
}

template <>
std::expected<std::pair<FuncF<double>, FuncFReference<double> >, Errors>
getFunctionFloat (const std::string_view &funcname)
{
  if (const auto it = find_function (funcF64, funcname))
    return std::make_pair (*it.value ()->func,
			   FuncFReference<double>{ *it.value ()->mpfrFunc });
  return std::unexpected (Errors::invalid_func);
}

template <>
std::expected<std::pair<FuncFF<float>, FuncFFReference<float> >, Errors>
getFunctionFloatFloat (const std::string_view &funcname)
{
  if (const auto it = find_function (funcF32F32, funcname))
    return std::make_pair (*it.value ()->func,
			   FuncFFReference<float>{ *it.value ()->mpfrFunc });
  return std::unexpected (Errors::invalid_func);
}

template <>
std::expected<std::pair<FuncFF<double>, FuncFFReference<double> >, Errors>
getFunctionFloatFloat (const std::string_view &funcname)
{
  if (const auto it = find_function (funcF64F64, funcname))
    return std::make_pair (*it.value ()->func,
			   FuncFFReference<double>{ *it.value ()->mpfrFunc });
  return std::unexpected (Errors::invalid_func);
}

template <>
std::expected<std::pair<FuncFpFp<float>, FuncFpFpReference<float> >, Errors>
getFunctionFloatpFloatp (const std::string_view &funcname)
{
  if (const auto it = find_function (funcFP32FP32, funcname))
    return std::make_pair (*it.value ()->func,
			   FuncFpFpReference<float>{ *it.value ()->mpfrFunc });
  return std::unexpected (Errors::invalid_func);
}

template <>
std::expected<std::pair<FuncFpFp<double>, FuncFpFpReference<double> >, Errors>
getFunctionFloatpFloatp (const std::string_view &funcname)
{
  if (const auto it = find_function (funcFP64FP64, funcname))
    return std::make_pair (*it.value ()->func, FuncFpFpReference<double>{
						   *it.value ()->mpfrFunc });
  return std::unexpected (Errors::invalid_func);
}

template <>
std::expected<std::pair<FuncFLLI<float>, FuncFLLIReference<float> >, Errors>
getFunctionFloatLLI (const std::string_view &funcname)
{
  if (const auto it = find_function (funcF32LLI, funcname))
    return std::make_pair (*it.value ()->func,
			   FuncFLLIReference<float>{ *it.value ()->mpfrFunc });
  return std::unexpected (Errors::invalid_func);
}

template <>
std::expected<std::pair<FuncFLLI<double>, FuncFLLIReference<double> >, Errors>
getFunctionFloatLLI (const std::string_view &funcname)
{
  if (const auto it = find_function (funcF64LLI, funcname))
    return std::make_pair (*it.value ()->func, FuncFLLIReference<double>{
						   *it.value ()->mpfrFunc });
  return std::unexpected (Errors::invalid_func);
}

std::expected<FunctionType, Errors>
getFunctionType (const std::string_view &funcname)
{
  if (contains_function (funcF32, funcname))
    return refimpls::FunctionType::f32_f;
  else if (contains_function (funcF32F32, funcname))
    return refimpls::FunctionType::f32_f_f;
  else if (contains_function (funcFP32FP32, funcname))
    return refimpls::FunctionType::f32_f_fp_fp;

  else if (contains_function (funcF64, funcname))
    return refimpls::FunctionType::f64_f;
  else if (contains_function (funcF64F64, funcname))
    return refimpls::FunctionType::f64_f_f;
  else if (contains_function (funcFP64FP64, funcname))
    return refimpls::FunctionType::f64_f_fp_fp;

  else if (contains_function (funcF32LLI, funcname))
    return refimpls::FunctionType::f32_f_lli;
  else if (contains_function (funcF64LLI, funcname))
    return refimpls::FunctionType::f64_f_lli;

  return std::unexpected (Errors::invalid_func);
}

} // namespace refimpls
