#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <functional>
#include <ranges>
#include <string>
#include <utility>

#include <fenv.h>
#include <mpfr.h>

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
init_ref_func<float> ()
{
  mpfr_set_emin (-148);
  mpfr_set_emax (128);
}

template <>
void
init_ref_func<double> ()
{
  mpfr_set_emin (-1073);
  mpfr_set_emax (1024);
}

extern "C"
{
#define _DEF_UNIVARIATE(__name)                                               \
  extern float cr_##__name##f (float);                                        \
  extern float ref_##__name##f (float, mpfr_rnd_t);                           \
  extern double cr_##__name (double);                                         \
  extern double ref_##__name (double, mpfr_rnd_t)

#define DEF_UNIVARIATE(__name) _DEF_UNIVARIATE (__name)

#define DEF_UNIVARIATE_WEAK(__name)                                           \
  extern float __name##f (float) __attribute__ ((weak));                      \
  extern double __name (double) __attribute__ ((weak));                       \
  _DEF_UNIVARIATE (__name)

#define _DEF_BIVARIATE(__name)                                                \
  extern float cr_##__name##f (float, float);                                 \
  extern float ref_##__name##f (float, float, mpfr_rnd_t);                    \
  extern double cr_##__name (double, double);                                 \
  extern double ref_##__name (double, double, mpfr_rnd_t)

#define DEF_BIVARIATE(__name) _DEF_BIVARIATE (__name)

#define DEF_BIVARIATE_WEAK(__name)                                            \
  extern float __name##f (float, float) __attribute__ ((weak, used));         \
  extern double __name (double, double) __attribute__ ((weak, used));         \
  _DEF_BIVARIATE (__name)

  DEF_BIVARIATE (atan2);
  DEF_UNIVARIATE_WEAK (atanpi);
  DEF_UNIVARIATE (acos);
  DEF_UNIVARIATE (acosh);
  DEF_UNIVARIATE_WEAK (acospi);
  DEF_UNIVARIATE (asin);
  DEF_UNIVARIATE (asinh);
  DEF_UNIVARIATE_WEAK (asinpi);
  DEF_UNIVARIATE (atan);
  DEF_UNIVARIATE (atanh);
  DEF_UNIVARIATE (cbrt);
  DEF_UNIVARIATE (cos);
  DEF_UNIVARIATE (cosh);
  DEF_UNIVARIATE_WEAK (cospi);
  DEF_UNIVARIATE (erf);
  DEF_UNIVARIATE (erfc);
  DEF_UNIVARIATE (exp);
  DEF_UNIVARIATE (expm1);
  DEF_UNIVARIATE_WEAK (exp10);
  DEF_UNIVARIATE_WEAK (exp10m1);
  DEF_UNIVARIATE (exp2);
  DEF_UNIVARIATE_WEAK (exp2m1);
  DEF_UNIVARIATE (lgamma);
  DEF_UNIVARIATE (log);
  DEF_UNIVARIATE (log1p);
  DEF_UNIVARIATE (log2);
  DEF_UNIVARIATE_WEAK (log2p1);
  DEF_UNIVARIATE (log10);
  DEF_UNIVARIATE_WEAK (log10p1);
  DEF_BIVARIATE (hypot);
  DEF_BIVARIATE (pow);
  DEF_UNIVARIATE_WEAK (rsqrt);
  DEF_UNIVARIATE (sin);
  DEF_UNIVARIATE (sinh);
  DEF_UNIVARIATE_WEAK (sinpi);
  DEF_UNIVARIATE (tan);
  DEF_UNIVARIATE (tanh);
  DEF_UNIVARIATE_WEAK (tanpi);
  DEF_UNIVARIATE (tgamma);

#undef _DEF_UNIVARIATE
#undef DEF_UNIVARIATE
#undef DEF_UNIVARIATE_WEAK
#undef _DEF_BIVARIATE
#undef DEF_BIVARIATE
#undef DEF_BIVARIATE_WEAK
};

typedef float (*binary32_univariate_t) (float);
typedef float (*binary32_univariate_mpfr_t) (float, mpfr_rnd_t);
typedef float (*binary32_bivariate_t) (float, float);
typedef float (*binary32_bivariate_mpfr_t) (float, float, mpfr_rnd_t);

typedef double (*binary64_univariate_t) (double);
typedef double (*binary64_univariate_mpfr_t) (double, mpfr_rnd_t);
typedef double (*binary64_bivariate_t) (double, double);
typedef double (*binary64_bivariate_mpfr_t) (double, double, mpfr_rnd_t);

template <typename T> struct ref_mode_univariate
{
  ref_mode_univariate (const T &func) : f (func) {}

  double
  operator() (double input, int rnd)
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

  T f;
};

template <typename T> struct ref_mode_bivariate
{
  ref_mode_bivariate (const T &func) : f (func) {}

  double
  operator() (double x, double y, int rnd)
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

  T f;
};

template <typename F, typename F_MPFR> struct univariate_functions_t
{
  const std::string name;
  volatile F func;
  F cr_func;
  F_MPFR mpfr_func;
};
typedef univariate_functions_t<binary32_univariate_t,
                               binary32_univariate_mpfr_t>
    binary32_univariate_functions_t;

static float
lgammaf_wrapper (float x)
{
  // the lgamma is not thread-safe
  int sign;
  return lgammaf_r (x, &sign);
}

// clang-format off
const static std::array binary32_univariate_functions = {
#define FUNC_DEF(name)                                                        \
  binary32_univariate_functions_t {                                           \
    #name, name, cr_##name, ref_##name                                        \
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
  binary32_univariate_functions_t {
    "lgammaf", lgammaf_wrapper, cr_lgammaf, ref_lgammaf
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

typedef univariate_functions_t<binary64_univariate_t,
                               binary64_univariate_mpfr_t>
    binary64_univariate_functions_t;

static double
lgamma_wrapper (double x)
{
  // the lgamma is not thread-safe
  int sign;
  return lgamma_r (x, &sign);
}

// clang-format off
const static std::array binary64_univariate_functions = {
#define FUNC_DEF(name)                                                        \
  binary64_univariate_functions_t {                                           \
    #name, name, cr_##name, ref_##name                                        \
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
  binary64_univariate_functions_t {
    "lgamma", lgamma_wrapper, cr_lgamma, ref_lgamma
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

template <typename F, typename F_MPFR> struct bivariate_functions_t
{
  const std::string name;
  F func;
  F cr_func;
  F_MPFR mpfr_func;
};

typedef univariate_functions_t<binary32_bivariate_t, binary32_bivariate_mpfr_t>
    binary32_bivariate_functions_t;

// clang-format off
const static std::array binary32_bivariate_functions = {
#define FUNC_DEF(name)                                                        \
  binary32_bivariate_functions_t {                                            \
    #name, name, cr_##name, ref_##name                                        \
  }
  FUNC_DEF (powf),
#undef FUNC_DEF
  };
// clang-format on

typedef univariate_functions_t<binary64_bivariate_t, binary64_bivariate_mpfr_t>
    binary64_bivariate_functions_t;

// clang-format off
const static std::array binary64_bivariate_functions = {
#define FUNC_DEF(name)                                                        \
  binary64_bivariate_functions_t {                                            \
    #name, name, cr_##name, ref_##name                                        \
  }
  FUNC_DEF (atan2),
  FUNC_DEF (hypot),
  FUNC_DEF (pow),
#undef FUNC_DEF
  };
// clang-format on

enum class function_error_t
{
  FUNCTION_NOT_FOUND
};

template <typename T, std::size_t N>
static std::expected<typename std::array<T, N>::const_iterator,
                     function_error_t>
find_function (const std::array<T, N> &funcs, const std::string_view &funcname)
{
  auto it = std::ranges::find_if (
      funcs, [&funcname] (const std::array<T, N>::const_reference &func) {
        return func.name == funcname;
      });
  if (it != funcs.end ())
    return it;
  return std::unexpected (function_error_t::FUNCTION_NOT_FOUND);
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

std::expected<std::pair<univariate_binary32_t, univariate_binary32_ref_t>,
              errors_t>
get_univariate_binary32 (const std::string_view &funcname, bool coremath)
{
  if (const auto it = find_function (binary32_univariate_functions, funcname))
    return std::make_pair (coremath ? *it.value ()->cr_func
                                    : *it.value ()->func,
                           ref_mode_univariate<binary32_univariate_mpfr_t>{
                               *it.value ()->mpfr_func });
  return std::unexpected (errors_t::invalid_func);
}

std::expected<std::pair<bivariate_binary32_t, bivariate_binary32_ref_t>,
              errors_t>
get_bivariate_binary32 (const std::string_view &funcname, bool coremath)
{
  if (const auto it = find_function (binary32_bivariate_functions, funcname))
    return std::make_pair (coremath ? *it.value ()->cr_func
                                    : *it.value ()->func,
                           ref_mode_bivariate<binary32_bivariate_mpfr_t>{
                               *it.value ()->mpfr_func });
  return std::unexpected (errors_t::invalid_func);
}

std::expected<std::pair<univariate_binary64_t, univariate_binary64_ref_t>,
              errors_t>
get_univariate_binary64 (const std::string_view &funcname, bool coremath)
{
  if (const auto it = find_function (binary64_univariate_functions, funcname))
    return std::make_pair (coremath ? *it.value ()->cr_func
                                    : *it.value ()->func,
                           ref_mode_univariate<binary64_univariate_mpfr_t>{
                               *it.value ()->mpfr_func });
  return std::unexpected (errors_t::invalid_func);
}

std::expected<std::pair<bivariate_binary64_t, bivariate_binary64_ref_t>,
              errors_t>
get_bivariate_binary64 (const std::string_view &funcname, bool coremath)
{
  if (const auto it = find_function (binary64_bivariate_functions, funcname))
    return std::make_pair (coremath ? *it.value ()->cr_func
                                    : *it.value ()->func,
                           ref_mode_bivariate<binary64_bivariate_mpfr_t>{
                               *it.value ()->mpfr_func });
  return std::unexpected (errors_t::invalid_func);
}

std::expected<func_type_t, errors_t>
get_func_type (const std::string_view &funcname)
{
  if (contains_function (binary32_univariate_functions, funcname))
    return refimpls::func_type_t::binary32_univariate;
  else if (contains_function (binary32_bivariate_functions, funcname))
    return refimpls::func_type_t::binary32_bivariate;
  if (contains_function (binary64_univariate_functions, funcname))
    return refimpls::func_type_t::binary64_univariate;
  else if (contains_function (binary64_bivariate_functions, funcname))
    return refimpls::func_type_t::binary64_bivariate;

  return std::unexpected (errors_t::invalid_func);
}

} // namespace refimpls
