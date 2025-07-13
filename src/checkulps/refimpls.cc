#include <algorithm>
#include <cmath>
#include <cstdint>
#include <string>
#include <utility>

#include <fenv.h>
#include <mpfr.h>

#include "refimpls.h"

namespace refimpls
{

extern "C"
{
#define _DEF_UNIVARIATE(__name)                                               \
  double cr_##__name (double);                                                \
  double ref_##__name (double, mpfr_rnd_t)

#define DEF_UNIVARIATE(__name) _DEF_UNIVARIATE (__name)

#define DEF_UNIVARIATE_WEAK(__name)                                           \
  double __name (double) __attribute__ ((weak));                              \
  _DEF_UNIVARIATE (__name)

#define _DEF_BIVARIATE(__name)                                                \
  double cr_##__name (double, double);                                        \
  double ref_##__name (double, double, mpfr_rnd_t)

#define DEF_BIVARIATE(__name) _DEF_BIVARIATE (__name)

#define DEF_BIVARIATE_WEAK(__name)                                            \
  double __name (double, double) __attribute__ ((weak, used));                \
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
  DEF_UNIVARIATE (exp10);
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

typedef double (*univariate_raw_t) (double);
typedef double (*univariate_mpfr_raw_t) (double, mpfr_rnd_t);
typedef double (*bivariate_raw_t) (double, double);
typedef double (*bivariate_mpfr_raw_t) (double, double, mpfr_rnd_t);

typedef std::function<double (double, mpfr_rnd_t)> univariate_mpfr_t;
typedef std::function<double (double, double, mpfr_rnd_t)> bivariate_mpfr_t;

struct ref_mode_univariate
{
  ref_mode_univariate (const univariate_mpfr_raw_t func) : f (func) {}

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

  const univariate_mpfr_t f;
};

struct ref_mode_bivariate
{
  ref_mode_bivariate (const bivariate_mpfr_raw_t &func) : f (func) {}

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

  const bivariate_mpfr_t f;
};

struct univariate_functions_t
{
  const std::string name;
  univariate_raw_t func;
  univariate_raw_t cr_func;
  univariate_mpfr_raw_t mpfr_func;
};

static double
lgamma_wrapper (double x)
{
  // the lgamma is not thread-safe
  int sign;
  return lgamma_r (x, &sign);
}

const static std::vector<univariate_functions_t> univariate_functions = {
#define FUNC_DEF(name)                                                        \
  {                                                                           \
    #name, name, cr_##name, ref_##name                                        \
  }
  FUNC_DEF (atanpi),  FUNC_DEF (acos),
  FUNC_DEF (acosh),   FUNC_DEF (acospi),
  FUNC_DEF (asin),    FUNC_DEF (asinh),
  FUNC_DEF (asinpi),  FUNC_DEF (atan),
  FUNC_DEF (atanh),   FUNC_DEF (cbrt),
  FUNC_DEF (cos),     FUNC_DEF (cosh),
  FUNC_DEF (cospi),   FUNC_DEF (erf),
  FUNC_DEF (erfc),    FUNC_DEF (exp),
  FUNC_DEF (expm1),   FUNC_DEF (exp10),
  FUNC_DEF (exp10m1), FUNC_DEF (exp2),
  FUNC_DEF (exp2m1),  { "lgamma", lgamma_wrapper, cr_lgamma, ref_lgamma },
  FUNC_DEF (log),     FUNC_DEF (log1p),
  FUNC_DEF (log2),    FUNC_DEF (log2p1),
  FUNC_DEF (log10),   FUNC_DEF (log10p1),
  FUNC_DEF (rsqrt),   FUNC_DEF (sin),
  FUNC_DEF (sinh),    FUNC_DEF (sinpi),
  FUNC_DEF (tan),     FUNC_DEF (tanh),
  FUNC_DEF (tanpi),   FUNC_DEF (tgamma),
#undef FUNC_DEF
};

struct bivariate_functions_t
{
  const std::string name;
  bivariate_raw_t func;
  bivariate_raw_t cr_func;
  bivariate_mpfr_raw_t mpfr_func;
};

const static std::vector<bivariate_functions_t> bivariate_functions = {
#define FUNC_DEF(name)                                                        \
  {                                                                           \
    #name, name, cr_##name, ref_##name                                        \
  }
  FUNC_DEF (atan2),
  FUNC_DEF (hypot),
#undef FUNC_DEF
};

static inline auto
find_univariate (const std::string_view &funcname)
{
  return std::find_if (univariate_functions.begin (),
                       univariate_functions.end (),
                       [&funcname] (const univariate_functions_t &func) {
                         return func.name == funcname;
                       });
}

static inline auto
find_bivariate (const std::string_view &funcname)
{
  return std::find_if (bivariate_functions.begin (),
                       bivariate_functions.end (),
                       [&funcname] (const bivariate_functions_t &func) {
                         return func.name == funcname;
                       });
}

std::expected<std::pair<univariate_t, univariate_ref_t>, errors_t>
get_univariate (const std::string_view &funcname, bool coremath)
{
  const auto it = find_univariate (funcname);
  if (it != univariate_functions.end ())
    return std::make_pair (coremath ? (*it).cr_func : (*it).func,
                           ref_mode_univariate{ (*it).mpfr_func });
  return std::unexpected (errors_t::invalid_func);
}

std::expected<std::pair<bivariate_t, bivariate_ref_t>, errors_t>
get_bivariate (const std::string_view &funcname, bool coremath)
{
  const auto it = find_bivariate (funcname);
  if (it != bivariate_functions.end ())
    return std::make_pair (coremath ? (*it).cr_func : (*it).func,
                           ref_mode_bivariate{ (*it).mpfr_func });
  return std::unexpected (errors_t::invalid_func);
}

std::expected<func_type_t, errors_t>
get_func_type (const std::string_view &funcname)
{
  if (find_univariate (funcname) != univariate_functions.end ())
    return refimpls::func_type_t::univariate;

  if (find_bivariate (funcname) != bivariate_functions.end ())
    return refimpls::func_type_t::bivariate;

  return std::unexpected (errors_t::invalid_func);
}

} // namespace refimpls
