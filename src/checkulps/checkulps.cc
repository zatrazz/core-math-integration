//
// Copyright (c) Adhemerval Zanella. All rights reserved.
// Licensed under the MIT license. See LICENSE file in the project root for
// details.
//

#include <algorithm>
#include <iostream>
#include <numbers>
#include <random>
#include <ranges>

#include <argparse/argparse.hpp>

#include <fenv.h>
#include <omp.h>

#include "description.h"
#include "floatranges.h"
#include "iohelper.h"
#include "refimpls.h"
#include "wyhash64.h"

// This is the threshold used by glibc that triggers a failure.
static constexpr auto k_max_ulp_str = "9.0";

using namespace refimpls;
using namespace iohelper;
typedef wyhash64 rng_t;

using clock_type = std::chrono::high_resolution_clock;

//
// issignaling: C11 macro that returns if a number is a signaling NaN.
//

template <std::floating_point T> bool issignaling (T);

template <>
bool
issignaling<float> (float x)
{
  union
  {
    float f;
    std::uint32_t u;
  } n = { .f = x };
  n.u ^= 0x00400000;
  return (n.u & 0x7fffffff) > 0x7fc00000;
}

template <>
bool
issignaling<double> (double x)
{
  union
  {
    double f;
    std::uint64_t u;
  } n = { .f = x };
  n.u ^= UINT64_C (0x0008000000000000);
  return (n.u & UINT64_C (0x7fffffffffffffff)) > UINT64_C (0x7ff8000000000000);
}

//
// round_mode_t: a class wrapper over C99 rounding modes, used to select which
//               one to test.  The default is to check for all rounding modes,
//               with the option to select a subset through command line.
//

struct round_mode_t
{
  const std::string name;
  const std::string abbrev;
  const int mode;

  constexpr
  round_mode_t (int rnd)
      : name (round_name (rnd)), abbrev (round_name (rnd, true)), mode (rnd)
  {
  }

  constexpr
  round_mode_t (const std::string &n, const std::string &a, int m)
      : name (n), abbrev (a), mode (m)
  {
  }

  bool
  operator== (const round_mode_t &other) const
  {
    return mode == other.mode;
  }

  bool
  operator== (const std::string &other) const
  {
    return abbrev == other;
  }

  bool
  operator== (int m) const
  {
    return mode == m;
  }

  constexpr static inline std::string
  round_name (int rnd, bool abbrev = false)
  {
    switch (rnd)
      {
      case FE_TONEAREST:
	return abbrev ? "rndn" : "FE_TONEAREST";
      case FE_UPWARD:
	return abbrev ? "rndu" : "FE_UPWARD";
      case FE_DOWNWARD:
	return abbrev ? "rndd" : "FE_DOWNWARD";
      case FE_TOWARDZERO:
	return abbrev ? "rndz" : "FE_TOWARDZERO";
      default:
	std::unreachable ();
      }
  }
};

typedef std::vector<round_mode_t> round_set;

// Use an array to keep insertion order, the idea is to first use the default
// rounding mode (FE_TONEAREST).
const static std::array k_round_modes
    = { round_mode_t (FE_TONEAREST), round_mode_t (FE_UPWARD),
	round_mode_t (FE_DOWNWARD), round_mode_t (FE_TOWARDZERO) };

const round_mode_t &
round_mode_from_rnd (int rnd)
{
  if (auto it = std::find (k_round_modes.begin (), k_round_modes.end (), rnd);
      it != k_round_modes.end ())
    return *it;
  std::abort ();
}

static constexpr std::vector<std::string>
split_with_ranges (const std::string_view &s, std::string_view delimiter)
{
  std::vector<std::string> tokens;
  for (const auto &subrange : s | std::views::split (delimiter))
    tokens.push_back (std::string (subrange.begin (), subrange.end ()));
  return tokens;
}

static const round_set
round_from_option (const std::string_view &rnds)
{
  round_set ret;

  auto round_modes = split_with_ranges (rnds, ",");
  for (auto &rnd : round_modes)
    if (auto it
	= std::find (k_round_modes.begin (), k_round_modes.end (), rnd);
	it != k_round_modes.end ())
      {
	if (std::ranges::contains (ret, *it))
	  error ("rounding mode already defined: {}", rnd);
	ret.push_back (*it);
      }
    else
      error ("invalid rounding mode: {}", rnd);

  return ret;
}

constexpr std::string
default_round_option ()
{
  return std::accumulate (
      std::next (k_round_modes.begin ()), k_round_modes.end (),
      k_round_modes.empty () ? ""
			     : std::string (k_round_modes.begin ()->abbrev),
      [] (const std::string &acc, const auto &rnd) {
	return acc + ',' + std::string (rnd.abbrev);
      });
}

//
// fail_mode_t: How to act when a failure is found:
//              - none:  report the found ULP distribution.
//              - first: print the error information and exit when first
//                       invalid or large error is found.
//              - all:   print the error information and continue checking.
//

enum class fail_mode_t
{
  none,
  first,
  all,
};

static const std::map<std::string_view, fail_mode_t> k_fail_modes
    = { { "none", fail_mode_t::none },
	{ "first", fail_mode_t::first },
	{ "all", fail_mode_t::all } };

fail_mode_t
fail_mode_from_options (const std::string_view &failmode)
{
  if (auto it = k_fail_modes.find (failmode); it != k_fail_modes.end ())
    return it->second;

  error ("invalid fail mode: {}", failmode);
}

// Returns the size of an ulp for VALUE.
template <typename F>
F
ulp (F value)
{
  F ulp;

  switch (std::fpclassify (value))
    {
    case FP_ZERO:
      /* Fall through...  */
    case FP_SUBNORMAL:
      ulp = std::ldexp (1.0, std::numeric_limits<F>::min_exponent
				 - std::numeric_limits<F>::digits);
      break;

    case FP_NORMAL:
      ulp = std::ldexp (1.0, std::ilogb (value)
				 - std::numeric_limits<F>::digits + 1);
      break;

    default:
      std::unreachable ();
      break;
    }
  return ulp;
}

/* Returns the number of ulps that GIVEN is away from EXPECTED.  */
template <typename F>
F
ulpdiff (F given, F expected)
{
  return std::fabs (given - expected) / ulp (expected);
}

template <typename F> using ulpacc_t = std::map<F, uint64_t>;

template <typename F>
static void
ulpacc_reduction (ulpacc_t<F> &inout, ulpacc_t<F> &in)
{
  for (const auto &ulp : in)
    inout[ulp.first] += ulp.second;
}

static int
get_max_thread (void)
{
#ifdef _OPENMP
  return omp_get_max_threads ();
#else
  return 1;
#endif
}

static int
get_thread_num (void)
{
#ifdef _OPENMP
  return omp_get_thread_num ();
#else
  return 0;
#endif
}

template <typename F>
static inline F
parse_range (const std::string &str)
{
  if (str == "-pi")
    return -std::numbers::pi_v<F>;
  else if (str == "pi")
    return std::numbers::pi_v<F>;
  else if (str == "2pi")
    return 2.0 * std::numbers::pi_v<F>;
  else if (str == "min")
    return std::numeric_limits<F>::min ();
  else if (str == "-min")
    return -std::numeric_limits<F>::min ();
  else if (str == "max")
    return std::numeric_limits<F>::max ();
  else if (str == "-max")
    return -std::numeric_limits<F>::max ();
  return std::stod (str);
}

//
// round_setup_t: Helper class to setup/reset the rounding mode, along with
//                extra setup required by MPFR.
//

template <typename F> class round_setup_t
{
  int saved_rnd;
  std::pair<mpfr_exp_t, mpfr_exp_t> saved_ref_param;

  void
  setup_rnd (int rnd)
  {
    if (fesetround (rnd) != 0)
      error ("fesetround ({}) failed (errno={})",
	     round_mode_t::round_name (rnd), errno);
  }

public:
  explicit round_setup_t (int rnd)
  {
    saved_rnd = fegetround ();
    setup_rnd (rnd);
    refimpls::setupReferenceImpl<F> ();
  }
  ~round_setup_t () { setup_rnd (saved_rnd); }

  // Prevent copying and moving to ensure single execution
  round_setup_t (const round_setup_t &) = delete;
  round_setup_t &operator= (const round_setup_t &) = delete;
  round_setup_t (round_setup_t &&) = delete;
  round_setup_t &operator= (round_setup_t &&) = delete;
};

// Accumulate histogram printer helpers.

template <typename F>
static void
print_acc (const std::string_view &rndname,
	   const Description::Sample1Arg<F> &sample, const ulpacc_t<F> &ulpacc)
{
  const std::uint64_t ulptotal = std::accumulate (
      ulpacc.begin (), ulpacc.end (), UINT64_C (0),
      [] (const uint64_t previous, const std::pair<double, uint64_t> &p) {
	return previous + p.second;
      });

  printlnTimestamp (
      "Checking rounding mode {:13}, range [{:9.2g},{:9.2g}], count {}",
      rndname, sample.arg.start, sample.arg.end, ulptotal);

  for (const auto &ulp : ulpacc)
    printlnTimestamp ("    {:g}: {:16} {:6.2f}%", ulp.first, ulp.second,
		      ((double) ulp.second / (double) ulptotal) * 100.0);
}

template <typename F>
static void
print_acc (const std::string_view &rndname,
	   const Description::Sample2Arg<F> &sample, const ulpacc_t<F> &ulpacc)
{
  const std::uint64_t ulptotal = std::accumulate (
      ulpacc.begin (), ulpacc.end (), UINT64_C (0),
      [] (const uint64_t previous, const std::pair<double, uint64_t> &p) {
	return previous + p.second;
      });

  printlnTimestamp ("Checking rounding mode {:13}, range x=[{:9.2g},{:9.2g}], "
		    "y=[{:9.2g},{:9.2g}], count {}",
		    rndname, sample.arg_x.start, sample.arg_x.end,
		    sample.arg_y.start, sample.arg_y.end, ulptotal);

  for (const auto &ulp : ulpacc)
    printlnTimestamp ("    {:g}: {:16} {:6.2f}%", ulp.first, ulp.second,
		      ((double) ulp.second / (double) ulptotal) * 100.0);
}

template <typename F>
static void
print_acc (const std::string_view &rndname,
	   const Description::Sample2ArgLli<F> &sample,
	   const ulpacc_t<F> &ulpacc)
{
  const std::uint64_t ulptotal = std::accumulate (
      ulpacc.begin (), ulpacc.end (), UINT64_C (0),
      [] (const uint64_t previous, const std::pair<double, uint64_t> &p) {
	return previous + p.second;
      });

  printlnTimestamp ("Checking rounding mode {:13}, range x=[{:9.2g},{:9.2g}], "
		    "y=[{},{}], count {}",
		    rndname, sample.arg_x.start, sample.arg_x.end,
		    sample.arg_y.start, sample.arg_y.end, ulptotal);

  for (const auto &ulp : ulpacc)
    printlnTimestamp ("    {:g}: {:16} {:6.2f}%", ulp.first, ulp.second,
		      ((double) ulp.second / (double) ulptotal) * 100.0);
}

template <typename F>
static void
print_acc (const std::string_view &rndname,
	   const Description::FullRange &sample, const ulpacc_t<F> &ulpacc)
{
  const std::uint64_t ulptotal = std::accumulate (
      ulpacc.begin (), ulpacc.end (), UINT64_C (0),
      [] (const uint64_t previous, const std::pair<double, uint64_t> &p) {
	return previous + p.second;
      });

  printlnTimestamp ("Checking rounding mode {:13}, {}", rndname, sample.name);

  for (const auto &ulp : ulpacc)
    printlnTimestamp ("    {:g}: {:16} {:6.2f}%", ulp.first, ulp.second,
		      ((double) ulp.second / (double) ulptotal) * 100.0);
}

static std::vector<rng_t::state_type> rng_states;

static void
init_random_state (void)
{
  std::random_device rd{ "/dev/random" };
  rng_states.resize (get_max_thread ());
  for (auto &s : rng_states)
    // std::random_device max is UINT32_MAX
    s = (rng_t::state_type) rd () << 32 | rd ();
}

template <typename F> struct result_t
{
  typedef F float_type;

  result_t (int r, F c, F e, F m)
      : rnd (round_mode_from_rnd (r)), computed (c), expected (e), max (m)
  {
    ulp = ulpdiff (computed, expected);
    if (std::isnan (ulp) || std::isinf (ulp))
      // Do not signal an error if the expected value is NaN/Inf.
      ulp = 0.0;

    // This avoid adding too much elements in the checking map for
    // implementation that has bad precision on non standard rounding
    // modes.
    if (ulp >= max)
      ulp = max;
  }

  bool
  check (void) const
  {
    return ulp < max;
  }

  bool
  check_full (void) const
  {
    if (issignaling (computed) || issignaling (expected))
      return false;
    else if (std::isnan (computed) && std::isnan (expected))
      return true;
    else if (std::isinf (computed) && std::isinf (expected))
      /* Test for sign of infinities.  */
      return std::signbit (computed) == std::signbit (expected);
    else if (std::isinf (computed) || std::isnan (computed)
	     || std::isinf (expected) || std::isnan (expected))
      return false;

    return check ();
  }

  virtual std::ostream &print (std::ostream &) const = 0;

  friend std::ostream &
  operator<< (std::ostream &os, const result_t &e)
  {
    return e.print (os);
  }

  const round_mode_t &rnd;
  F computed;
  F expected;
  F ulp;
  F max;
};

template <typename F> struct result_f_t : public result_t<F>
{
  explicit result_f_t (int r, F i, F c, F e, F m)
      : result_t<F> (r, c, e, m), input (i)
  {
  }

  virtual std::ostream &
  print (std::ostream &os) const
  {
    os << std::format (
	"{} ulp={:1.0f} input=0x{:a} computed=0x{:a} expected=0x{:a}\n",
	result_t<F>::rnd.name, result_t<F>::ulp, input, result_t<F>::computed,
	result_t<F>::expected);
    return os;
  }

  F input;
};

template <typename F> struct result_f_fp_fp_t
{
  typedef F float_type;

  explicit result_f_fp_fp_t (int r, F i, F c1, F c2, F e1, F e2, F m)
      : rnd (round_mode_from_rnd (r)), input (i), computed1 (c1),
	computed2 (c2), expected1 (e1), expected2 (e2), max (m)
  {
    F ulp0 = calc_ulp (computed1, expected1);
    F ulp1 = calc_ulp (computed2, expected2);
    ulp = ulp1 > ulp0 ? ulp1 : ulp0;
  }

  inline F
  calc_ulp (F computed, F expected)
  {
    F ulp = ulpdiff (computed, expected);
    if (std::isnan (ulp) || std::isinf (ulp))
      // Do not signal an error if the expected value is NaN/Inf.
      return 0.0;
    return ulp >= max ? max : ulp;
  }

  bool
  check (void) const
  {
    return ulp < max;
  }

  bool
  check_full (void) const
  {
    if (issignaling (computed1) || issignaling (computed2)
	|| issignaling (expected2) || issignaling (expected2))
      return false;

    else if ((std::isnan (computed1) && std::isnan (expected1))
	     && (std::isnan (computed2) && std::isnan (expected2)))
      return true;

    else if ((std::isinf (computed1) && std::isinf (expected1))
	     && (std::isinf (computed2) && std::isinf (expected2)))
      /* Test for sign of infinities.  */
      return std::signbit (computed1) == std::signbit (expected1)
	     && std::signbit (computed2) == std::signbit (expected2);

    return check ();
  }

  friend std::ostream &
  operator<< (std::ostream &os, const result_f_fp_fp_t &e)
  {
    return e.print (os);
  }

  virtual std::ostream &
  print (std::ostream &os) const
  {
    os << std::format ("{} ulp={:1.0f} input=0x{:a} computed=(0x{:a} 0x{:a}) "
		       "expected=(0x{:a} 0x{:a}i)\n",
		       rnd.name, ulp, input, computed1, computed2, expected1,
		       expected2);
    return os;
  }

  const round_mode_t &rnd;
  F input;
  F computed1;
  F computed2;
  F expected1;
  F expected2;
  F ulp;
  F max;
};

template <typename F> struct result_f_f_t : public result_t<F>
{
  explicit result_f_f_t (int r, F i0, F i1, F c, F e, F m)
      : result_t<F> (r, c, e, m), input0 (i0), input1 (i1)
  {
  }

  virtual std::ostream &
  print (std::ostream &os) const
  {
    os << std::format ("{} ulp={:1.0f} input=(0x{:a},0x{:a}) computed=0x{:a} "
		       "expected=0x{:a}\n",
		       result_t<F>::rnd.name, result_t<F>::ulp, input0, input1,
		       result_t<F>::computed, result_t<F>::expected);
    return os;
  }

  F input0;
  F input1;
};

template <typename F> struct result_f_lli_t : public result_t<F>
{
public:
  explicit result_f_lli_t (int r, F i0, long long int i1, F c, F e, F m)
      : result_t<F> (r, c, e, m), input0 (i0), input1 (i1)
  {
  }

  virtual std::ostream &
  print (std::ostream &os) const
  {
    os << std::format ("{} ulp={:1.0f} input=(0x{:a},{}) computed=0x{:a} "
		       "expected=0x{:a}\n",
		       result_t<F>::rnd.name, result_t<F>::ulp, input0, input1,
		       result_t<F>::computed, result_t<F>::expected);
    return os;
  }

  F input0;
  long long int input1;
};

template <typename RET> struct sample_random_base_f_t
{
  virtual std::unique_ptr<RET>
  operator() (rng_t &,
	      std::uniform_real_distribution<typename RET::float_type> &,
	      int) const
      = 0;
};

template <typename FUNC, typename FUNC_REF, typename RET>
class sample_random_f_t : public sample_random_base_f_t<RET>
{
  typedef typename RET::float_type float_type;

  const FUNC &func;
  const FUNC_REF &ref_func;
  float_type max_ulp;

public:
  sample_random_f_t (FUNC &f, FUNC_REF &ref_f, RET::float_type mulp)
      : func (f), ref_func (ref_f), max_ulp (mulp)
  {
  }

  std::unique_ptr<RET>
  operator() (rng_t &gen, std::uniform_real_distribution<float_type> &dist,
	      int rnd) const
  {
    float_type input = dist (gen);
    float_type computed = func (input);
    float_type expected = ref_func (input, rnd);

    return std::make_unique<RET> (rnd, input, computed, expected, max_ulp);
  }
};
template <typename F>
using random_f_t
    = sample_random_f_t<FuncF<F>, FuncFReference<F>, result_f_t<F> >;

template <typename FUNC, typename FUNC_REF, typename RET>
class sample_random_f_fp_fp_t : public sample_random_base_f_t<RET>
{
  typedef typename RET::float_type float_type;

  const FUNC &func;
  const FUNC_REF &ref_func;
  float_type max_ulp;

public:
  sample_random_f_fp_fp_t (FUNC &f, FUNC_REF &ref_f, RET::float_type mulp)
      : func (f), ref_func (ref_f), max_ulp (mulp)
  {
  }

  std::unique_ptr<RET>
  operator() (rng_t &gen, std::uniform_real_distribution<float_type> &dist,
	      int rnd) const
  {
    float_type input = dist (gen);

    float_type computed0, computed1;
    func (input, &computed0, &computed1);

    float_type expected0, expected1;
    ref_func (input, &expected0, &expected1, rnd);

    return std::make_unique<RET> (rnd, input, computed0, computed1, expected0,
				  expected1, max_ulp);
  }
};
template <typename F>
using random_f_fp_fp_t
    = sample_random_f_fp_fp_t<FuncFpFp<F>, FuncFpFpReference<F>,
			      result_f_fp_fp_t<F> >;

template <typename RET> struct sample_random_base_f_f_t
{
  virtual std::unique_ptr<RET> operator() (
      rng_t &, std::uniform_real_distribution<typename RET::float_type> &,
      std::uniform_real_distribution<typename RET::float_type> &, int) const
      = 0;
};

template <typename FUNC, typename FUNC_REF, typename RET>
class sample_random_f_f_t : public sample_random_base_f_f_t<RET>
{
  typedef typename RET::float_type float_type;

  const FUNC &func;
  const FUNC_REF &ref_func;
  float_type max_ulp;

public:
  sample_random_f_f_t (FUNC &f, FUNC_REF &ref_f, RET::float_type mulp)
      : func (f), ref_func (ref_f), max_ulp (mulp)
  {
  }

  std::unique_ptr<RET>
  operator() (rng_t &gen, std::uniform_real_distribution<float_type> &dist_x,
	      std::uniform_real_distribution<float_type> &dist_y,
	      int rnd) const
  {
    float_type input0 = dist_x (gen);
    float_type input1 = dist_y (gen);
    float_type computed = func (input0, input1);
    float_type expected = ref_func (input0, input1, rnd);

    return std::make_unique<RET> (rnd, input0, input1, computed, expected,
				  max_ulp);
  }
};
template <typename F>
using random_f_f_t
    = sample_random_f_f_t<FuncFF<F>, FuncFFReference<F>, result_f_f_t<F> >;

template <typename RET> struct sample_random_base_f_lli_t
{
  virtual std::unique_ptr<RET>
  operator() (rng_t &,
	      std::uniform_real_distribution<typename RET::float_type> &,
	      std::uniform_int_distribution<long long int> &, int) const
      = 0;
};

template <typename FUNC, typename FUNC_REF, typename RET>
class sample_random_f_lli_t : public sample_random_base_f_lli_t<RET>
{
  typedef typename RET::float_type float_type;

  const FUNC &func;
  const FUNC_REF &ref_func;
  float_type max_ulp;

public:
  sample_random_f_lli_t (FUNC &f, FUNC_REF &ref_f, RET::float_type mulp)
      : func (f), ref_func (ref_f), max_ulp (mulp)
  {
  }

  std::unique_ptr<RET>
  operator() (rng_t &gen, std::uniform_real_distribution<float_type> &distx,
	      std::uniform_int_distribution<long long int> &disty,
	      int rnd) const
  {
    float_type input0 = distx (gen);
    long long int input1 = disty (gen);
    float_type computed = func (input0, input1);
    float_type expected = ref_func (input0, input1, rnd);

    return std::make_unique<RET> (rnd, input0, input1, computed, expected,
				  max_ulp);
  }
};
template <typename F>
using random_f_lli_t = sample_random_f_lli_t<FuncFLLI<F>, FuncFLLIReference<F>,
					     result_f_lli_t<F> >;

template <typename RET> struct sample_full_t
{
  RET::float_type max_ulp;

  sample_full_t (RET::float_type maxp) : max_ulp (maxp) {}

  virtual std::unique_ptr<RET> operator() (uint64_t, int) const = 0;
};

template <typename FUNC, typename FUNC_REF, typename RET>
class sample_full_f_t : public sample_full_t<RET>
{
  typedef typename RET::float_type float_type;

  const FUNC &func;
  const FUNC_REF &ref_func;

public:
  sample_full_f_t (FUNC &f, FUNC_REF &ref_f, RET::float_type mulp)
      : sample_full_t<RET> (mulp), func (f), ref_func (ref_f)
  {
  }

  std::unique_ptr<RET>
  operator() (uint64_t n, int rnd) const
  {
    float_type input = floatrange::Limits<float_type>::from (n);
    float_type computed = func (input);
    float_type expected = ref_func (input, rnd);

    return std::make_unique<RET> (rnd, input, computed, expected,
				  sample_full_t<RET>::max_ulp);
  }
};
template <typename F>
using full_f_t = sample_full_f_t<FuncF<F>, FuncFReference<F>, result_f_t<F> >;

template <typename FUNC, typename FUNC_REF, typename RET>
class sample_full_f_fp_fp_t : public sample_full_t<RET>
{
  typedef typename RET::float_type float_type;

  const FUNC &func;
  const FUNC_REF &ref_func;

public:
  sample_full_f_fp_fp_t (FUNC &f, FUNC_REF &ref_f, RET::float_type mulp)
      : sample_full_t<RET> (mulp), func (f), ref_func (ref_f)
  {
  }

  std::unique_ptr<RET>
  operator() (uint64_t n, int rnd) const
  {
    float_type input = floatrange::Limits<float_type>::from (n);
    float_type computed0, computed1;
    func (input, &computed0, &computed1);

    float_type expected0, expected1;
    ref_func (input, &expected0, &expected1, rnd);

    return std::make_unique<RET> (rnd, input, computed0, computed1, expected0,
				  expected1, sample_full_t<RET>::max_ulp);
  }
};
template <typename F>
using full_f_fp_fp_t = sample_full_f_fp_fp_t<FuncFpFp<F>, FuncFpFpReference<F>,
					     result_f_fp_fp_t<F> >;

template <typename RET>
static void
check_random_f (
    const std::string_view &funcname, const sample_random_base_f_t<RET> &funcs,
    const Description::Sample1Arg<typename RET::float_type> &sample,
    const round_set &round_modes, fail_mode_t failmode)
{
  using float_type = typename RET::float_type;

  std::vector<rng_t> gens (rng_states.size ());

  for (auto &rnd : round_modes)
    {
      /* Reseed the RNG generator with the same seed to check the same numbers
	 for each rounding mode.  */
      for (unsigned i = 0; i < rng_states.size (); i++)
	gens[i] = rng_t (rng_states[i]);

#pragma omp declare reduction(                                                \
	ulpacc_reduction : ulpacc_t<float_type> : ulpacc_reduction(           \
		omp_out, omp_in))                                             \
    initializer(omp_priv = ulpacc_t<float_type> ())

      auto start = clock_type::now ();

      std::uniform_real_distribution<float_type> dist (sample.arg.start,
						       sample.arg.end);

      ulpacc_t<float_type> ulpaccrange;

#pragma omp parallel firstprivate(dist, failmode) shared(sample, rnd)
      {
	round_setup_t<float_type> round_setup (rnd.mode);

#pragma omp for reduction(ulpacc_reduction : ulpaccrange)
	for (std::uint64_t i = 0; i < sample.count; i++)
	  {
	    auto ret = funcs (gens[get_thread_num ()], dist, rnd.mode);
	    if (!ret->check ())
	      switch (failmode)
		{
		case fail_mode_t::first:
#pragma omp critical
		  std::cerr << *ret;
		  std::exit (EXIT_FAILURE);
		  break;
		case fail_mode_t::all:
#pragma omp critical
		  std::cerr << *ret;
		  break;
		}
	    ulpaccrange[ret->ulp] += 1;
	  }
      }

      print_acc (rnd.name, sample, ulpaccrange);

      auto end = clock_type::now ();
      printlnTimestamp (
	  "Elapsed time {}",
	  std::chrono::duration_cast<std::chrono::duration<double> > (
	      end - start));
      printlnTimestamp ("");
    }
}

template <typename RET>
static void
check_random_f_f (
    const std::string_view &funcname,
    const sample_random_base_f_f_t<RET> &funcs,
    const Description::Sample2Arg<typename RET::float_type> &sample,
    const round_set &round_modes, fail_mode_t failmode)
{
  using float_type = typename RET::float_type;

  std::vector<rng_t> gens (rng_states.size ());

  for (auto &rnd : round_modes)
    {
      /* Reseed the RNG generator with the same seed to check the same numbers
	 for each rounding mode.  */
      for (unsigned i = 0; i < rng_states.size (); i++)
	gens[i] = rng_t (rng_states[i]);

#pragma omp declare reduction(                                                \
	ulpacc_reduction : ulpacc_t<float_type> : ulpacc_reduction(           \
		omp_out, omp_in))                                             \
    initializer(omp_priv = ulpacc_t<float_type> ())

      auto start = clock_type::now ();

      std::uniform_real_distribution<float_type> dist_x (sample.arg_x.start,
							 sample.arg_x.end);
      std::uniform_real_distribution<float_type> dist_y (sample.arg_y.start,
							 sample.arg_y.end);

      ulpacc_t<float_type> ulpaccrange;

#pragma omp parallel firstprivate(dist_x, dist_y, failmode) shared(sample, rnd)
      {
	round_setup_t<float_type> round_setup (rnd.mode);

#pragma omp for reduction(ulpacc_reduction : ulpaccrange)
	for (std::uint64_t i = 0; i < sample.count; i++)
	  {
	    auto ret
		= funcs (gens[get_thread_num ()], dist_x, dist_y, rnd.mode);
	    if (!ret->check ())
	      switch (failmode)
		{
		case fail_mode_t::first:
#pragma omp critical
		  std::cerr << *ret;
		  std::exit (EXIT_FAILURE);
		  break;
		case fail_mode_t::all:
#pragma omp critical
		  std::cerr << *ret;
		  break;
		}
	    ulpaccrange[ret->ulp] += 1;
	  }
      }

      print_acc (rnd.name, sample, ulpaccrange);

      auto end = clock_type::now ();
      printlnTimestamp (
	  "Elapsed time {}",
	  std::chrono::duration_cast<std::chrono::duration<double> > (
	      end - start));
      printlnTimestamp ("");
    }
}

template <typename RET>
static void
check_random_f_lli (
    const std::string_view &funcname,
    const sample_random_base_f_lli_t<RET> &funcs,
    const Description::Sample2ArgLli<typename RET::float_type> &sample,
    const round_set &round_modes, fail_mode_t failmode)
{
  using float_type = typename RET::float_type;
  using y_type = long long int;

  std::vector<rng_t> gens (rng_states.size ());

  for (auto &rnd : round_modes)
    {
      /* Reseed the RNG generator with the same seed to check the same numbers
	 for each rounding mode.  */
      for (unsigned i = 0; i < rng_states.size (); i++)
	gens[i] = rng_t (rng_states[i]);

#pragma omp declare reduction(                                                \
	ulpacc_reduction : ulpacc_t<float_type> : ulpacc_reduction(           \
		omp_out, omp_in))                                             \
    initializer(omp_priv = ulpacc_t<float_type> ())

      auto start = clock_type::now ();

      std::uniform_real_distribution<float_type> dist_x (sample.arg_x.start,
							 sample.arg_x.end);
      std::uniform_int_distribution<y_type> dist_y (sample.arg_y.start,
						    sample.arg_y.end);

      ulpacc_t<float_type> ulpaccrange;

#pragma omp parallel firstprivate(dist_x, dist_y, failmode) shared(sample, rnd)
      {
	round_setup_t<float_type> round_setup (rnd.mode);

#pragma omp for reduction(ulpacc_reduction : ulpaccrange)
	for (std::uint64_t i = 0; i < sample.count; i++)
	  {
	    auto ret
		= funcs (gens[get_thread_num ()], dist_x, dist_y, rnd.mode);
	    if (!ret->check ())
	      switch (failmode)
		{
		case fail_mode_t::first:
#pragma omp critical
		  std::cerr << *ret;
		  std::exit (EXIT_FAILURE);
		  break;
		case fail_mode_t::all:
#pragma omp critical
		  std::cerr << *ret;
		  break;
		}
	    ulpaccrange[ret->ulp] += 1;
	  }
      }

      print_acc (rnd.name, sample, ulpaccrange);

      auto end = clock_type::now ();
      printlnTimestamp (
	  "Elapsed time {}",
	  std::chrono::duration_cast<std::chrono::duration<double> > (
	      end - start));
      printlnTimestamp ("");
    }
}

template <typename RET>
static void
check_full_f (const std::string_view &funcname,
	      const sample_full_t<RET> &funcs,
	      const Description::FullRange &sample,
	      const round_set &round_modes, fail_mode_t failmode)
{
  using float_type = typename RET::float_type;

  for (auto &rnd : round_modes)
    {
#pragma omp declare reduction(                                                \
	ulpacc_reduction : ulpacc_t<float_type> : ulpacc_reduction(           \
		omp_out, omp_in))                                             \
    initializer(omp_priv = ulpacc_t<float_type> ())

      ulpacc_t<float_type> ulpaccrange;

#pragma omp parallel firstprivate(failmode) shared(funcs, rnd)
      {
	round_setup_t<float_type> round_setup (rnd.mode);

// Out of range inputs might take way less time than normal one, also use
// a large chunk size to minimize the overhead from dynamic scheduline.
#pragma omp for reduction(ulpacc_reduction : ulpaccrange) schedule(dynamic)
	for (std::uint64_t i = sample.start; i < sample.end; i++)
	  {
	    auto ret = funcs (i, rnd.mode);
	    if (!ret->check_full ())
	      switch (failmode)
		{
		case fail_mode_t::first:
#pragma omp critical
		  std::cerr << *ret;
		  std::exit (EXIT_FAILURE);
		  break;
		case fail_mode_t::all:
#pragma omp critical
		  std::cerr << *ret;
		  break;
		}
	    ulpaccrange[ret->ulp] += 1;
	  }
      }

      print_acc (rnd.name, sample, ulpaccrange);
      printlnTimestamp ("");
    }
}

template <typename F>
static void
run_f (const Description &desc, const round_set &round_modes,
       fail_mode_t failmode, const std::string &max_ulp_str)
{
  auto func = get_f<F> (desc.FunctionName).value ();
  if (!func.first)
    error ("libc does not provide {}", desc.FunctionName);

  const auto max_ulp = floatrange::fromStr<F> (max_ulp_str);
  if (!max_ulp)
    error ("invalid floating point: {}", max_ulp_str);

  printlnTimestamp ("Checking function {}", desc.FunctionName);
  printlnTimestamp ("");

  auto start = clock_type::now ();

  for (auto &sample : desc.Samples)
    {
      if (auto *psample = std::get_if<Description::Sample1Arg<F> > (&sample))
	check_random_f (
	    desc.FunctionName,
	    random_f_t<F>{ func.first, func.second, max_ulp.value () },
	    *psample, round_modes, failmode);
      else if (auto *psample = std::get_if<Description::FullRange> (&sample))
	check_full_f (desc.FunctionName,
		      full_f_t<F>{ func.first, func.second, max_ulp.value () },
		      *psample, round_modes, failmode);
      else
	error ("invalid sample type");
    }

  auto end = clock_type::now ();
  printlnTimestamp (
      "Total elapsed time {}",
      std::chrono::duration_cast<std::chrono::duration<double> > (end
								  - start));
}

template <typename F>
static void
run_f_fp_fp (const Description &desc, const round_set &round_modes,
	     fail_mode_t failmode, const std::string &max_ulp_str)
{
  auto func = get_f_fp_fp<F> (desc.FunctionName).value ();
  if (!func.first)
    error ("libc does not provide {}", desc.FunctionName);

  const auto max_ulp = floatrange::fromStr<F> (max_ulp_str);
  if (!max_ulp)
    error ("invalid floating point: {}", max_ulp_str);

  printlnTimestamp ("Checking function {}", desc.FunctionName);
  printlnTimestamp ("");

  auto start = clock_type::now ();

  for (auto &sample : desc.Samples)
    {
      if (auto *psample = std::get_if<Description::Sample1Arg<F> > (&sample))
	check_random_f (
	    desc.FunctionName,
	    random_f_fp_fp_t<F>{ func.first, func.second, max_ulp.value () },
	    *psample, round_modes, failmode);
      else if (auto *psample = std::get_if<Description::FullRange> (&sample))
	check_full_f (
	    desc.FunctionName,
	    full_f_fp_fp_t<F>{ func.first, func.second, max_ulp.value () },
	    *psample, round_modes, failmode);
      else
	error ("invalid sample type");
    }

  auto end = clock_type::now ();
  printlnTimestamp (
      "Total elapsed time {}",
      std::chrono::duration_cast<std::chrono::duration<double> > (end
								  - start));
}

template <typename F>
static void
run_f_f (const Description &desc, const round_set &round_modes,
	 fail_mode_t failmode, const std::string &max_ulp_str)
{
  auto func = get_f_f<F> (desc.FunctionName).value ();
  if (!func.first)
    error ("libc does not provide {}", desc.FunctionName);

  const auto max_ulp = floatrange::fromStr<F> (max_ulp_str);
  if (!max_ulp)
    error ("invalid floating point: {}", max_ulp_str);

  printlnTimestamp ("Checking FunctionName {}", desc.FunctionName);
  printlnTimestamp ("");

  auto start = clock_type::now ();

  for (auto &sample : desc.Samples)
    {
      if (auto *psample = std::get_if<Description::Sample2Arg<F> > (&sample))
	check_random_f_f (
	    desc.FunctionName,
	    random_f_f_t<F>{ func.first, func.second, max_ulp.value () },
	    *psample, round_modes, failmode);
      else
	error ("invalid sample type");
    }

  auto end = clock_type::now ();
  printlnTimestamp (
      "Total elapsed time {}",
      std::chrono::duration_cast<std::chrono::duration<double> > (end
								  - start));
}

template <typename F>
static void
run_f_lli (const Description &desc, const round_set &round_modes,
	   fail_mode_t failmode, const std::string &max_ulp_str)
{
  auto func = get_f_lli<F> (desc.FunctionName).value ();
  if (!func.first)
    error ("libc does not provide {}", desc.FunctionName);

  const auto max_ulp = floatrange::fromStr<F> (max_ulp_str);
  if (!max_ulp)
    error ("invalid floating point: {}", max_ulp_str);

  printlnTimestamp ("Checking FunctionName {}", desc.FunctionName);
  printlnTimestamp ("");

  auto start = clock_type::now ();

  for (auto &sample : desc.Samples)
    {
      if (auto *psample
	  = std::get_if<Description::Sample2ArgLli<F> > (&sample))
	check_random_f_lli (
	    desc.FunctionName,
	    random_f_lli_t<F>{ func.first, func.second, max_ulp.value () },
	    *psample, round_modes, failmode);
      else
	error ("invalid sample type");
    }

  auto end = clock_type::now ();
  printlnTimestamp (
      "Total elapsed time {}",
      std::chrono::duration_cast<std::chrono::duration<double> > (end
								  - start));
}

int
main (int argc, char *argv[])
{
  argparse::ArgumentParser options ("checkulps");

  options.add_argument ("--description", "-d")
      .help ("input JSON descriptiorn file")
      .required ();

  options.add_argument ("--rounding", "-r")
      .help ("rounding modes to test")
      .default_value (default_round_option ());

  options.add_argument ("--failure", "-f")
      .help ("failure mode")
      .default_value ("none");

  options.add_argument ("--maxulps", "-m")
      .help ("max ULP used in check")
      .default_value (k_max_ulp_str);

  try
    {
      options.parse_args (argc, argv);
    }
  catch (const std::runtime_error &err)
    {
      error (std::string (err.what ()));
    }

  const round_set round_modes
      = round_from_option (options.get<std::string> ("-r"));

  fail_mode_t failmode
      = fail_mode_from_options (options.get<std::string> ("-f"));

  Description desc;
  if (auto r = desc.parse (options.get<std::string> ("-d")); !r)
    error ("{}", r.error ());

  std::string max_ulp = options.get<std::string> ("-m");

  init_random_state ();

  auto functype = get_func_type (desc.FunctionName);
  if (!functype)
    error ("invalid FunctionName: {}", desc.FunctionName);

  switch (functype.value ())
    {
    case refimpls::func_type_t::f32_f:
      run_f<float> (desc, round_modes, failmode, max_ulp);
      break;
    case refimpls::func_type_t::f64_f:
      run_f<double> (desc, round_modes, failmode, max_ulp);
      break;

    case refimpls::func_type_t::f32_f_f:
      run_f_f<float> (desc, round_modes, failmode, max_ulp);
      break;
    case refimpls::func_type_t::f64_f_f:
      run_f_f<double> (desc, round_modes, failmode, max_ulp);
      break;

    case refimpls::func_type_t::f32_f_lli:
      run_f_lli<float> (desc, round_modes, failmode, max_ulp);
      break;
    case refimpls::func_type_t::f64_f_lli:
      run_f_lli<double> (desc, round_modes, failmode, max_ulp);
      break;

    case refimpls::func_type_t::f32_f_fp_fp:
      run_f_fp_fp<float> (desc, round_modes, failmode, max_ulp);
      break;
    case refimpls::func_type_t::f64_f_fp_fp:
      run_f_fp_fp<double> (desc, round_modes, failmode, max_ulp);
      break;

    default:
      error ("function type \"{}\" not implemented", functype.value ());
    }
}
