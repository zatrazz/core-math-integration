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

  double acospi (double) __attribute__ ((weak)); 
  double asinpi (double) __attribute__ ((weak)); 
};;

typedef std::function<double(double, mpfr_rnd_t rnd)> univariate_mpfr_t;

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

struct funcs_t
{
  const std::string name;
  univariate_t      func;
  univariate_t      cr_func;
  univariate_mpfr_t mpfr_func;
};

const static std::vector<funcs_t> math_functions = {
#define FUNC_DEF(name) { #name, name, cr_ ## name, ref_ ## name }
  FUNC_DEF (asin),
  FUNC_DEF (asinh),
  FUNC_DEF (asinpi),

  FUNC_DEF (acos),
  FUNC_DEF (acosh),
  FUNC_DEF (acospi),
#undef FUNC_DEF
};

std::expected<univariate_t, errors_t>
get_univariate (const std::string_view& funcname, bool coremath)
{
  const auto it = std::find_if(math_functions.begin(), math_functions.end(),
			       [&funcname](const funcs_t &func) {
			         return func.name == funcname;
				});
  if (it != math_functions.end())
    return coremath ? (*it).cr_func : (*it).func;
  return std::unexpected (errors_t::invalid_func);
}

std::expected<univariate_ref_t, errors_t>
get_univariate_ref (const std::string_view &funcname)
{
  const auto it = std::find_if(math_functions.begin(), math_functions.end(),
			       [&funcname](const funcs_t &func) {
			         return func.name == funcname;
				});
  if (it != math_functions.end())
    return ref_mode_univariate {(*it).mpfr_func};
  return std::unexpected (errors_t::invalid_func);
}

}
