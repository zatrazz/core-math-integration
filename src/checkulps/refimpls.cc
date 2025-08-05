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
setup_ref_impl<float> ()
{
  mpfr_set_emin (-148);
  mpfr_set_emax (128);
}

template <>
void
setup_ref_impl<double> ()
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

template <typename F, typename F_MPFR> struct func_f_desc_t
{
  const std::string name;
  volatile F func;
  F_MPFR mpfr_func;
};
typedef func_f_desc_t<func_f_t<float>, func_f_mpfr_t<float> >
    func_f32_f_desc_t;

static float
lgammaf_wrapper (float x)
{
  // The lgamma is not thread-safe, and we can ignore the sign.
  int sign;
  return lgammaf_r (x, &sign);
}

// clang-format off
const static std::array func_f32_f_desc = {
#define FUNC_DEF(name)                                                        \
  func_f32_f_desc_t {                                           \
    #name, name, ref_##name                                        \
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
  func_f32_f_desc_t {
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

typedef func_f_desc_t<func_f_t<double>, func_f_mpfr_t<double> >
    func_f64_f_desc_t;

static double
lgamma_wrapper (double x)
{
  // The lgamma is not thread-safe, and we can ignore the sign.
  int sign;
  return lgamma_r (x, &sign);
}

// clang-format off
const static std::array func_f64_f_desc = {
#define FUNC_DEF(name)                                                        \
  func_f64_f_desc_t {                                           \
    #name, name, ref_##name                                        \
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
  func_f64_f_desc_t {
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

template <typename F, typename F_MPFR> struct func_f_f_desc_t
{
  const std::string name;
  F func;
  F_MPFR mpfr_func;
};

typedef func_f_desc_t<func_f_f_t<float>, func_f_f_mpfr_t<float> >
    func_f32_f_f_desc_t;

// clang-format off
const static std::array func_f32_f_f_desc = {
#define FUNC_DEF(name)                                                        \
  func_f32_f_f_desc_t {                                            \
    #name, name, ref_##name                                        \
  }
  FUNC_DEF (atan2f),
  FUNC_DEF (hypotf),
  FUNC_DEF (powf),
  FUNC_DEF (powrf),
#undef FUNC_DEF
  };
// clang-format on

typedef func_f_desc_t<func_f_f_t<double>, func_f_f_mpfr_t<double> >
    func_f64_f_f_desc_t;

// clang-format off
const static std::array func_f64_f_f_desc = {
#define FUNC_DEF(name)                                                        \
  func_f64_f_f_desc_t {                                            \
    #name, name, ref_##name                                        \
  }
  FUNC_DEF (atan2),
  FUNC_DEF (hypot),
  FUNC_DEF (pow),
  FUNC_DEF (powr),
#undef FUNC_DEF
  };
// clang-format on

template <typename F, typename F_MPFR> struct func_f_lli_desc_t
{
  const std::string name;
  F func;
  F_MPFR mpfr_func;
};

typedef func_f_desc_t<func_f_lli_t<float>, func_f_lli_mpfr_t<float> >
    func_f32_f_lli_desc_t;

// clang-format off
const static std::array func_f32_f_lli_desc = {
#define FUNC_DEF(name)                                                        \
  func_f32_f_lli_desc_t {                                            \
    #name, name, ref_##name                                        \
  }
  FUNC_DEF (compoundnf),
  FUNC_DEF (pownf),
  FUNC_DEF (rootnf),
#undef FUNC_DEF
  };
// clang-format on

typedef func_f_desc_t<func_f_lli_t<double>, func_f_lli_mpfr_t<double> >
    func_f64_f_lli_desc_t;

// clang-format off
const static std::array func_f64_f_lli_desc = {
#define FUNC_DEF(name)                                                        \
  func_f64_f_lli_desc_t {                                            \
    #name, name, ref_##name                                        \
  }
  FUNC_DEF (compoundn),
  FUNC_DEF (pown),
  FUNC_DEF (rootn),
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
  auto it = std::find_if (
      funcs.begin (), funcs.end (),
      [&funcname] (const std::array<T, N>::const_reference &func) {
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

template <>
std::expected<std::pair<func_f_t<float>, func_ref_t<float> >, errors_t>
get_f (const std::string_view &funcname)
{
  if (const auto it = find_function (func_f32_f_desc, funcname))
    return std::make_pair (*it.value ()->func,
                           func_ref_t<float>{ *it.value ()->mpfr_func });
  return std::unexpected (errors_t::invalid_func);
}

template <>
std::expected<std::pair<func_f_t<double>, func_ref_t<double> >, errors_t>
get_f (const std::string_view &funcname)
{
  if (const auto it = find_function (func_f64_f_desc, funcname))
    return std::make_pair (*it.value ()->func,
                           func_ref_t<double>{ *it.value ()->mpfr_func });
  return std::unexpected (errors_t::invalid_func);
}

template <>
std::expected<std::pair<func_f_f_t<float>, func_f_f_ref_t<float> >, errors_t>
get_f_f (const std::string_view &funcname)
{
  if (const auto it = find_function (func_f32_f_f_desc, funcname))
    return std::make_pair (*it.value ()->func,
                           func_f_f_ref_t<float>{ *it.value ()->mpfr_func });
  return std::unexpected (errors_t::invalid_func);
}

template <>
std::expected<std::pair<func_f_f_t<double>, func_f_f_ref_t<double> >, errors_t>
get_f_f (const std::string_view &funcname)
{
  if (const auto it = find_function (func_f64_f_f_desc, funcname))
    return std::make_pair (*it.value ()->func,
                           func_f_f_ref_t<double>{ *it.value ()->mpfr_func });
  return std::unexpected (errors_t::invalid_func);
}

template <>
std::expected<std::pair<func_f_lli_t<float>, func_f_lli_ref_t<float> >,
              errors_t>
get_f_lli (const std::string_view &funcname)
{
  if (const auto it = find_function (func_f32_f_lli_desc, funcname))
    return std::make_pair (*it.value ()->func,
                           func_f_lli_ref_t<float>{ *it.value ()->mpfr_func });
  return std::unexpected (errors_t::invalid_func);
}

template <>
std::expected<std::pair<func_f_lli_t<double>, func_f_lli_ref_t<double> >,
              errors_t>
get_f_lli (const std::string_view &funcname)
{
  if (const auto it = find_function (func_f64_f_lli_desc, funcname))
    return std::make_pair (*it.value ()->func, func_f_lli_ref_t<double>{
                                                   *it.value ()->mpfr_func });
  return std::unexpected (errors_t::invalid_func);
}

std::expected<func_type_t, errors_t>
get_func_type (const std::string_view &funcname)
{
  if (contains_function (func_f32_f_desc, funcname))
    return refimpls::func_type_t::f32_f;
  else if (contains_function (func_f32_f_f_desc, funcname))
    return refimpls::func_type_t::f32_f_f;

  else if (contains_function (func_f64_f_desc, funcname))
    return refimpls::func_type_t::f64_f;
  else if (contains_function (func_f64_f_f_desc, funcname))
    return refimpls::func_type_t::f64_f_f;

  else if (contains_function (func_f32_f_lli_desc, funcname))
    return refimpls::func_type_t::f32_f_lli;
  else if (contains_function (func_f64_f_lli_desc, funcname))
    return refimpls::func_type_t::f64_f_lli;

  return std::unexpected (errors_t::invalid_func);
}

} // namespace refimpls
