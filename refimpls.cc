#include <cmath>
#include <cstdint>
#include <string>
#include <utility>
#include <algorithm>

#include <fenv.h>
#include <mpfr.h>

#include "refimpls.h"

namespace refimpls
{

extern "C" {
  // From core-math
  double cr_asin (double);
  double cr_asinpi (double);
  double cr_asinh (double);
  double cr_acos (double);
  double cr_acospi (double);
  double cr_acosh (double);
  double cr_atan (double);
  double cr_atanh (double);
  double cr_atanpi (double);
  double cr_cos (double);
  double cr_sin (double);

  double acospi (double) __attribute__ ((weak)); 
  double asinpi (double) __attribute__ ((weak)); 
  double atanpi (double) __attribute__ ((weak));
};;

typedef std::function<double(double, mpfr_rnd_t rnd)> univariate_mpfr_t;
typedef std::function<double(double, double, mpfr_rnd_t rnd)> bivariate_mpfr_t;

struct ref_mode_univariate
{
  ref_mode_univariate(const univariate_mpfr_t& func) : f(func) {}

  double operator()(double input, int rnd)
  {
    switch (rnd)
    {
    case FE_TONEAREST:  return f (input, MPFR_RNDN);
    case FE_UPWARD:     return f (input, MPFR_RNDU);
    case FE_DOWNWARD:   return f (input, MPFR_RNDD);
    case FE_TOWARDZERO: return f (input, MPFR_RNDZ);
    default:	        std::unreachable();
    };
  }

  const univariate_mpfr_t &f;
};

struct ref_mode_bivariate
{
  ref_mode_bivariate(const bivariate_mpfr_t& func) : f(func) {}

  double operator()(double x, double y, int rnd)
  {
    switch (rnd)
    {
    case FE_TONEAREST:  return f (x, y, MPFR_RNDN);
    case FE_UPWARD:     return f (x, y, MPFR_RNDU);
    case FE_DOWNWARD:   return f (x, y, MPFR_RNDD);
    case FE_TOWARDZERO: return f (x, y, MPFR_RNDZ);
    default:	        std::unreachable();
    };
  }

  const bivariate_mpfr_t &f;
};

static double
ref_acos (double x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_init2 (y, 53);
  mpfr_set_d (y, x, MPFR_RNDN);
  int inex = mpfr_acos (y, y, rnd);
  mpfr_subnormalize (y, inex, rnd);
  double ret = mpfr_get_d (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

/* code from MPFR */
double
ref_acospi (double x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_init2 (y, 53);
  mpfr_set_d (y, x, MPFR_RNDN);
  mpfr_acospi (y, y, rnd);
  /* no need to call mpfr_subnormalize since the smallest non-zero value
     is obtained for x=0x1.fffffffffffffp-1, and is 0x1.45f306dc9c883p-28
     (rounded to nearest) */
  double ret = mpfr_get_d (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

typedef union {double f; std::uint64_t u;} b64u64_u;

double
ref_acosh (double x, mpfr_rnd_t rnd)
{
  b64u64_u ix = {.f = x};
  if(__builtin_expect(ix.u<=0x3ff0000000000000ull, 0)){
    if(ix.u==0x3ff0000000000000ull) return 0;
    return __builtin_nan("x<1");
  }
  if(__builtin_expect(ix.u>=0x7ff0000000000000ull, 0)){
    uint64_t aix = ix.u<<1;
    if(ix.u==0x7ff0000000000000ull || aix>((uint64_t)0x7ff<<53)) return x; // +inf or nan
    return __builtin_nan("x<1");
  }

  mpfr_t y;
  mpfr_init2 (y, 53);
  mpfr_set_d (y, x, MPFR_RNDN);
  mpfr_acosh (y, y, rnd);
  double ret = mpfr_get_d(y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

double
ref_asin (double x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_init2 (y, 53);
  mpfr_set_d (y, x, MPFR_RNDN);
  int inex = mpfr_asin (y, y, rnd);
  mpfr_subnormalize (y, inex, rnd);
  double ret = mpfr_get_d (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

double
ref_asinpi (double x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_init2 (y, 53);
  mpfr_set_d (y, x, MPFR_RNDN);
  int inex = mpfr_asinpi (y, y, rnd);
  mpfr_subnormalize (y, inex, rnd);
  double ret = mpfr_get_d (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

double ref_asinh (double x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_init2(y, 53);
  mpfr_set_d(y, x, MPFR_RNDN);
  int d = mpfr_asinh(y, y, rnd);
  mpfr_subnormalize(y, d, rnd);
  double ret = mpfr_get_d(y, MPFR_RNDN);
  mpfr_clear(y);
  return ret;
}

double ref_atan (double x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_init2 (y, 53);
  mpfr_set_d (y, x, MPFR_RNDN);
  int inex = mpfr_atan (y, y, rnd);
  mpfr_subnormalize (y, inex, rnd);
  double ret = mpfr_get_d (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

double ref_atanh (double x,  mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_init2(y, 53);
  mpfr_set_d(y, x, MPFR_RNDN);
  int inex = mpfr_atanh(y, y, rnd);
  mpfr_subnormalize(y, inex, rnd);
  double r = mpfr_get_d(y, MPFR_RNDN);
  mpfr_clear(y);
  return r;
}

double
ref_atanpi (double x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_init2 (y, 53);
  mpfr_set_d (y, x, MPFR_RNDN);
  int inex = mpfr_atanpi (y, y, rnd);
  mpfr_subnormalize (y, inex, rnd);
  double ret = mpfr_get_d (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

double ref_sin (double x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_init2 (y, 53);
  mpfr_set_d (y, x, MPFR_RNDN);
  int inex = mpfr_sin (y, y, rnd);
  mpfr_subnormalize (y, inex, rnd);
  double ret = mpfr_get_d (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

double ref_cos (double x, mpfr_rnd_t rnd)
{
  mpfr_t y;
  mpfr_init2 (y, 53);
  mpfr_set_d (y, x, MPFR_RNDN);
  int inex = mpfr_cos (y, y, rnd);
  mpfr_subnormalize (y, inex, rnd);
  double ret = mpfr_get_d (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}

struct univariate_functions_t
{
  const std::string name;
  univariate_t      func;
  univariate_t      cr_func;
  univariate_mpfr_t mpfr_func;
};

const static std::vector<univariate_functions_t> univariate_functions = {
#define FUNC_DEF(name) { #name, name, cr_ ## name, ref_ ## name }
  FUNC_DEF (asin),
  FUNC_DEF (asinh),
  FUNC_DEF (asinpi),

  FUNC_DEF (acos),
  FUNC_DEF (acosh),
  FUNC_DEF (acospi),

  FUNC_DEF (atan),
  FUNC_DEF (atanh),
  FUNC_DEF (atanpi),

  FUNC_DEF (cos),
  FUNC_DEF (sin),
#undef FUNC_DEF
};

struct bivariate_functions_t
{
  const std::string name;
  bivariate_t      func;
  bivariate_t      cr_func;
  bivariate_mpfr_t mpfr_func;
};

const static std::vector<bivariate_functions_t> bivariate_functions = {
#define FUNC_DEF(name) { #name, name, cr_ ## name, ref_ ## name }
#undef FUNC_DEF
};

static inline auto find_univariate (const std::string_view& funcname)
{
  return std::find_if (univariate_functions.begin(),
		       univariate_functions.end(),
		       [&funcname](const univariate_functions_t& func) {
		         return func.name == funcname;
		       });
}

static inline auto find_bivariate (const std::string_view& funcname)
{
  return std::find_if (bivariate_functions.begin(),
		       bivariate_functions.end(),
		       [&funcname](const bivariate_functions_t& func) {
		         return func.name == funcname;
		       });
}

std::expected<std::pair<univariate_t, univariate_ref_t>, errors_t>
get_univariate (const std::string_view& funcname, bool coremath)
{
  const auto it = find_univariate (funcname);
  if (it != univariate_functions.end())
    return std::make_pair(coremath ? (*it).cr_func : (*it).func,
			  ref_mode_univariate {(*it).mpfr_func});
  return std::unexpected (errors_t::invalid_func);
}

std::expected<std::pair<bivariate_t, bivariate_ref_t>, errors_t>
get_bivariate (const std::string_view& funcname, bool coremath)
{
  const auto it = find_bivariate (funcname);
  if (it != bivariate_functions.end())
    return std::make_pair(coremath ? (*it).cr_func : (*it).func,
			  ref_mode_bivariate {(*it).mpfr_func});
  return std::unexpected (errors_t::invalid_func);
}

std::expected<func_type_t, errors_t>
get_func_type (const std::string_view& funcname)
{
  if (find_univariate (funcname) != univariate_functions.end())
    return refimpls::func_type_t::univariate;

  if (find_bivariate (funcname) != bivariate_functions.end())
    return refimpls::func_type_t::bivariate;

  return std::unexpected (errors_t::invalid_func);
}

} //namespace refimpls
