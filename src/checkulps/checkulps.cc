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
#include "strhelper.h"

// This is the threshold used by glibc that triggers a failure.
static constexpr auto k_max_ulp_str = "9.0";

using namespace refimpls;
using namespace iohelper;
typedef wyhash64 RngType;

using ClockType = std::chrono::high_resolution_clock;

//
// isSignaling: C11 macro that returns if a number is a signaling NaN.
//

template <std::floating_point T> bool isSignaling (T);

template <>
bool
isSignaling<float> (float x)
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
isSignaling<double> (double x)
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
// RoundMode: a class wrapper over C99 rounding modes, used to select which
//               one to test.  The default is to check for all rounding modes,
//               with the option to select a subset through command line.
//

struct RoundMode
{
  const std::string name;
  const std::string abbrev;
  const int mode;

  constexpr
  RoundMode (int m)
      : name (roundName (m)), abbrev (roundName (m, true)), mode (m)
  {
  }

  constexpr
  RoundMode (const std::string &n, const std::string &a, int m)
      : name (n), abbrev (a), mode (m)
  {
  }

  bool
  operator== (const RoundMode &other) const
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
  roundName (int mode, bool abbrev = false)
  {
    switch (mode)
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

typedef std::vector<RoundMode> RoundSet;

// Use an array to keep insertion order, the idea is to first use the default
// rounding mode (FE_TONEAREST).
static const std::array kRoundModes
    = { RoundMode (FE_TONEAREST), RoundMode (FE_UPWARD),
	RoundMode (FE_DOWNWARD), RoundMode (FE_TOWARDZERO) };

const RoundMode &
roundModeFromRound (int rnd)
{
  if (auto it = std::find (kRoundModes.begin (), kRoundModes.end (), rnd);
      it != kRoundModes.end ())
    return *it;
  std::abort ();
}

static const RoundSet
roundFromOption (const std::string_view &rnds)
{
  RoundSet ret;

  auto roundModes = strhelper::splitWithRanges (rnds, ",");
  for (auto &rnd : roundModes)
    if (auto it = std::find (kRoundModes.begin (), kRoundModes.end (), rnd);
	it != kRoundModes.end ())
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
defaultRoundOption ()
{
  return std::accumulate (
      std::next (kRoundModes.begin ()), kRoundModes.end (),
      kRoundModes.empty () ? "" : std::string (kRoundModes.begin ()->abbrev),
      [] (const std::string &acc, const auto &rnd) {
	return acc + ',' + std::string (rnd.abbrev);
      });
}

//
// FailMode: How to act when a failure is found:
//              - none:  report the found ULP distribution.
//              - first: print the error information and exit when first
//                       invalid or large error is found.
//              - all:   print the error information and continue checking.
//

enum class FailMode
{
  NONE,
  FIRST,
  ALL,
};

static const std::map<std::string_view, FailMode> failModes
    = { { "none", FailMode::NONE },
	{ "first", FailMode::FIRST },
	{ "all", FailMode::ALL } };

FailMode
failModeFromOptions (const std::string_view &failmode)
{
  if (auto it = failModes.find (failmode); it != failModes.end ())
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

template <typename F> using UlpAccumulator = std::map<F, uint64_t>;

template <typename F>
static void
ulpAccumulatorReduction (UlpAccumulator<F> &inout, UlpAccumulator<F> &in)
{
  for (const auto &ulp : in)
    inout[ulp.first] += ulp.second;
}

static int
getMaxThread (void)
{
#ifdef _OPENMP
  return omp_get_max_threads ();
#else
  return 1;
#endif
}

static int
getThreadNum (void)
{
#ifdef _OPENMP
  return omp_get_thread_num ();
#else
  return 0;
#endif
}

//
// RoundSetup: Helper class to setup/reset the rounding mode, along with
//                extra setup required by MPFR.
//

template <typename F> class RoundSetup
{
  int savedMode;
  std::pair<mpfr_exp_t, mpfr_exp_t> savedRefParam;

  void
  setup (int mode)
  {
    if (fesetround (mode) != 0)
      error ("fesetround ({}) failed (errno={})", RoundMode::roundName (mode),
	     errno);
  }

public:
  explicit RoundSetup (int mode)
  {
    savedMode = fegetround ();
    setup (mode);
    refimpls::setupReferenceImpl<F> ();
  }
  ~RoundSetup () { setup (savedMode); }

  // Prevent copying and moving to ensure single execution
  RoundSetup (const RoundSetup &) = delete;
  RoundSetup &operator= (const RoundSetup &) = delete;
  RoundSetup (RoundSetup &&) = delete;
  RoundSetup &operator= (RoundSetup &&) = delete;
};

// Accumulate histogram printer helpers.

template <typename F>
static void
printAccumulator (const std::string_view &rndname,
		  const Description::Sample1Arg<F> &sample,
		  const UlpAccumulator<F> &ulpacc)
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
printAccumulator (const std::string_view &rndname,
		  const Description::Sample2Arg<F> &sample,
		  const UlpAccumulator<F> &ulpacc)
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
printAccumulator (const std::string_view &rndname,
		  const Description::Sample2ArgLli<F> &sample,
		  const UlpAccumulator<F> &ulpacc)
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
printAccumulator (const std::string_view &rndname,
		  const Description::FullRange &sample,
		  const UlpAccumulator<F> &ulpacc)
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

static std::vector<RngType::state_type> rngStates;

static void
initRandomState (void)
{
  std::random_device rd{ "/dev/random" };
  rngStates.resize (getMaxThread ());
  for (auto &s : rngStates)
    // std::random_device max is UINT32_MAX
    s = (RngType::state_type) rd () << 32 | rd ();
}

template <typename F> struct Result
{
  typedef F FloatType;

  Result (int r, F c, F e, F m)
      : roundMode (roundModeFromRound (r)), computed (c), expected (e), max (m)
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
  checkFull (void) const
  {
    if (isSignaling (computed) || issignaling (expected))
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
  operator<< (std::ostream &os, const Result &e)
  {
    return e.print (os);
  }

  const RoundMode &roundMode;
  FloatType computed;
  FloatType expected;
  FloatType ulp;
  FloatType max;
};

template <typename F> struct ResultFloat : public Result<F>
{
  explicit ResultFloat (int r, F i, F c, F e, F m)
      : Result<F> (r, c, e, m), input (i)
  {
  }

  virtual std::ostream &
  print (std::ostream &os) const
  {
    os << std::format (
	"{} ulp={:1.0f} input=0x{:a} computed=0x{:a} expected=0x{:a}\n",
	Result<F>::roundMode.name, Result<F>::ulp, input, Result<F>::computed,
	Result<F>::expected);
    return os;
  }

  F input;
};

template <typename F> struct ResultFloatpFloatp
{
  typedef F FloatType;

  explicit ResultFloatpFloatp (int r, F i, F c1, F c2, F e1, F e2, F m)
      : roundMode (roundModeFromRound (r)), input (i), computed1 (c1),
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
  checkFull (void) const
  {
    if (isSignaling (computed1) || issignaling (computed2)
	|| isSignaling (expected2) || issignaling (expected2))
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
  operator<< (std::ostream &os, const ResultFloatpFloatp &e)
  {
    return e.print (os);
  }

  virtual std::ostream &
  print (std::ostream &os) const
  {
    os << std::format ("{} ulp={:1.0f} input=0x{:a} computed=(0x{:a} 0x{:a}) "
		       "expected=(0x{:a} 0x{:a}i)\n",
		       roundMode.name, ulp, input, computed1, computed2,
		       expected1, expected2);
    return os;
  }

  const RoundMode &roundMode;
  FloatType input;
  FloatType computed1;
  FloatType computed2;
  FloatType expected1;
  FloatType expected2;
  FloatType ulp;
  FloatType max;
};

template <typename F> struct ResultFloatFloat : public Result<F>
{
  explicit ResultFloatFloat (int r, F i0, F i1, F c, F e, F m)
      : Result<F> (r, c, e, m), input0 (i0), input1 (i1)
  {
  }

  virtual std::ostream &
  print (std::ostream &os) const
  {
    os << std::format ("{} ulp={:1.0f} input=(0x{:a},0x{:a}) computed=0x{:a} "
		       "expected=0x{:a}\n",
		       Result<F>::roundMode.name, Result<F>::ulp, input0,
		       input1, Result<F>::computed, Result<F>::expected);
    return os;
  }

  F input0;
  F input1;
};

template <typename F> struct ResultFloatLLI : public Result<F>
{
public:
  explicit ResultFloatLLI (int r, F i0, long long int i1, F c, F e, F m)
      : Result<F> (r, c, e, m), input0 (i0), input1 (i1)
  {
  }

  virtual std::ostream &
  print (std::ostream &os) const
  {
    os << std::format ("{} ulp={:1.0f} input=(0x{:a},{}) computed=0x{:a} "
		       "expected=0x{:a}\n",
		       Result<F>::roundMode.name, Result<F>::ulp, input0,
		       input1, Result<F>::computed, Result<F>::expected);
    return os;
  }

  F input0;
  long long int input1;
};

template <typename RET> struct SampleRandomBase
{
  virtual std::unique_ptr<RET>
  operator() (RngType &,
	      std::uniform_real_distribution<typename RET::FloatType> &,
	      int) const
      = 0;
};

template <typename FUNC, typename FUNC_REF, typename RET>
class SampleRandomFloat : public SampleRandomBase<RET>
{
  typedef typename RET::FloatType FloatType;

  const FUNC &func;
  const FUNC_REF &ref_func;
  FloatType max_ulp;

public:
  SampleRandomFloat (FUNC &f, FUNC_REF &ref_f, RET::FloatType mulp)
      : func (f), ref_func (ref_f), max_ulp (mulp)
  {
  }

  std::unique_ptr<RET>
  operator() (RngType &gen, std::uniform_real_distribution<FloatType> &dist,
	      int rnd) const
  {
    FloatType input = dist (gen);
    FloatType computed = func (input);
    FloatType expected = ref_func (input, rnd);

    return std::make_unique<RET> (rnd, input, computed, expected, max_ulp);
  }
};
template <typename F>
using RandomFloat
    = SampleRandomFloat<FuncF<F>, FuncFReference<F>, ResultFloat<F> >;

template <typename FUNC, typename FUNC_REF, typename RET>
class SampleRandomFloatpFloatp : public SampleRandomBase<RET>
{
  typedef typename RET::FloatType FloatType;

  const FUNC &func;
  const FUNC_REF &ref_func;
  FloatType max_ulp;

public:
  SampleRandomFloatpFloatp (FUNC &f, FUNC_REF &ref_f, RET::FloatType mulp)
      : func (f), ref_func (ref_f), max_ulp (mulp)
  {
  }

  std::unique_ptr<RET>
  operator() (RngType &gen, std::uniform_real_distribution<FloatType> &dist,
	      int rnd) const
  {
    FloatType input = dist (gen);

    FloatType computed0, computed1;
    func (input, &computed0, &computed1);

    FloatType expected0, expected1;
    ref_func (input, &expected0, &expected1, rnd);

    return std::make_unique<RET> (rnd, input, computed0, computed1, expected0,
				  expected1, max_ulp);
  }
};
template <typename F>
using RandomFloatpFloatp
    = SampleRandomFloatpFloatp<FuncFpFp<F>, FuncFpFpReference<F>,
			       ResultFloatpFloatp<F> >;

template <typename RET> struct SampleRandomFloatFloatBase
{
  virtual std::unique_ptr<RET> operator() (
      RngType &, std::uniform_real_distribution<typename RET::FloatType> &,
      std::uniform_real_distribution<typename RET::FloatType> &, int) const
      = 0;
};

template <typename FUNC, typename FUNC_REF, typename RET>
class SampleRandomFloatFloat : public SampleRandomFloatFloatBase<RET>
{
  typedef typename RET::FloatType FloatType;

  const FUNC &func;
  const FUNC_REF &ref_func;
  FloatType max_ulp;

public:
  SampleRandomFloatFloat (FUNC &f, FUNC_REF &ref_f, RET::FloatType mulp)
      : func (f), ref_func (ref_f), max_ulp (mulp)
  {
  }

  std::unique_ptr<RET>
  operator() (RngType &gen, std::uniform_real_distribution<FloatType> &distX,
	      std::uniform_real_distribution<FloatType> &distY, int rnd) const
  {
    FloatType input0 = distX (gen);
    FloatType input1 = distY (gen);
    FloatType computed = func (input0, input1);
    FloatType expected = ref_func (input0, input1, rnd);

    return std::make_unique<RET> (rnd, input0, input1, computed, expected,
				  max_ulp);
  }
};
template <typename F>
using RandomFloatFloat = SampleRandomFloatFloat<FuncFF<F>, FuncFFReference<F>,
						ResultFloatFloat<F> >;

template <typename RET> struct SampleRandomFloatLLIBase
{
  virtual std::unique_ptr<RET>
  operator() (RngType &,
	      std::uniform_real_distribution<typename RET::FloatType> &,
	      std::uniform_int_distribution<long long int> &, int) const
      = 0;
};

template <typename FUNC, typename FUNC_REF, typename RET>
class SampleRandomFloatLLI : public SampleRandomFloatLLIBase<RET>
{
  typedef typename RET::FloatType FloatType;

  const FUNC &func;
  const FUNC_REF &ref_func;
  FloatType max_ulp;

public:
  SampleRandomFloatLLI (FUNC &f, FUNC_REF &ref_f, RET::FloatType mulp)
      : func (f), ref_func (ref_f), max_ulp (mulp)
  {
  }

  std::unique_ptr<RET>
  operator() (RngType &gen, std::uniform_real_distribution<FloatType> &distx,
	      std::uniform_int_distribution<long long int> &disty,
	      int rnd) const
  {
    FloatType input0 = distx (gen);
    long long int input1 = disty (gen);
    FloatType computed = func (input0, input1);
    FloatType expected = ref_func (input0, input1, rnd);

    return std::make_unique<RET> (rnd, input0, input1, computed, expected,
				  max_ulp);
  }
};
template <typename F>
using RandomFloatLLI = SampleRandomFloatLLI<FuncFLLI<F>, FuncFLLIReference<F>,
					    ResultFloatLLI<F> >;

template <typename RET> struct SampleFull
{
  RET::FloatType max_ulp;

  SampleFull (RET::FloatType maxp) : max_ulp (maxp) {}

  virtual std::unique_ptr<RET> operator() (uint64_t, int) const = 0;
};

template <typename FUNC, typename FUNC_REF, typename RET>
class SampleFullFloat : public SampleFull<RET>
{
  typedef typename RET::FloatType FloatType;

  const FUNC &func;
  const FUNC_REF &ref_func;

public:
  SampleFullFloat (FUNC &f, FUNC_REF &ref_f, RET::FloatType mulp)
      : SampleFull<RET> (mulp), func (f), ref_func (ref_f)
  {
  }

  std::unique_ptr<RET>
  operator() (uint64_t n, int rnd) const
  {
    FloatType input = floatrange::Limits<FloatType>::from (n);
    FloatType computed = func (input);
    FloatType expected = ref_func (input, rnd);

    return std::make_unique<RET> (rnd, input, computed, expected,
				  SampleFull<RET>::max_ulp);
  }
};
template <typename F>
using FullFloat
    = SampleFullFloat<FuncF<F>, FuncFReference<F>, ResultFloat<F> >;

template <typename FUNC, typename FUNC_REF, typename RET>
class SampleFullFloatpFloatp : public SampleFull<RET>
{
  typedef typename RET::FloatType FloatType;

  const FUNC &func;
  const FUNC_REF &ref_func;

public:
  SampleFullFloatpFloatp (FUNC &f, FUNC_REF &ref_f, RET::FloatType mulp)
      : SampleFull<RET> (mulp), func (f), ref_func (ref_f)
  {
  }

  std::unique_ptr<RET>
  operator() (uint64_t n, int rnd) const
  {
    FloatType input = floatrange::Limits<FloatType>::from (n);
    FloatType computed0, computed1;
    func (input, &computed0, &computed1);

    FloatType expected0, expected1;
    ref_func (input, &expected0, &expected1, rnd);

    return std::make_unique<RET> (rnd, input, computed0, computed1, expected0,
				  expected1, SampleFull<RET>::max_ulp);
  }
};
template <typename F>
using FullFloatpFloatp
    = SampleFullFloatpFloatp<FuncFpFp<F>, FuncFpFpReference<F>,
			     ResultFloatpFloatp<F> >;

template <typename RET>
static void
checkRandomFloat (
    const std::string_view &funcname, const SampleRandomBase<RET> &funcs,
    const Description::Sample1Arg<typename RET::FloatType> &sample,
    const RoundSet &roundModes, FailMode failmode)
{
  using FloatType = typename RET::FloatType;

  std::vector<RngType> gens (rngStates.size ());

  for (auto &rnd : roundModes)
    {
      /* Reseed the RNG generator with the same seed to check the same numbers
	 for each rounding mode.  */
      for (unsigned i = 0; i < rngStates.size (); i++)
	gens[i] = RngType (rngStates[i]);

#pragma omp declare reduction(                                                \
	ulpAccumulatorReduction : UlpAccumulator<                             \
		FloatType> : ulpAccumulatorReduction(omp_out, omp_in))        \
    initializer(omp_priv = UlpAccumulator<FloatType> ())

      auto start = ClockType::now ();

      std::uniform_real_distribution<FloatType> dist (sample.arg.start,
						      sample.arg.end);

      UlpAccumulator<FloatType> ulpaccrange;

#pragma omp parallel firstprivate(dist, failmode) shared(sample, rnd)
      {
	RoundSetup<FloatType> roundSetup (rnd.mode);

#pragma omp for reduction(ulpAccumulatorReduction : ulpaccrange)
	for (std::uint64_t i = 0; i < sample.count; i++)
	  {
	    auto ret = funcs (gens[getThreadNum ()], dist, rnd.mode);
	    if (!ret->check ())
	      switch (failmode)
		{
		case FailMode::FIRST:
#pragma omp critical
		  std::cerr << *ret;
		  std::exit (EXIT_FAILURE);
		  break;
		case FailMode::ALL:
#pragma omp critical
		  std::cerr << *ret;
		  [[fallthrough]];
		default:
		  break;
		}
	    ulpaccrange[ret->ulp] += 1;
	  }
      }

      printAccumulator (rnd.name, sample, ulpaccrange);

      auto end = ClockType::now ();
      printlnTimestamp (
	  "Elapsed time {}",
	  std::chrono::duration_cast<std::chrono::duration<double> > (
	      end - start));
      printlnTimestamp ("");
    }
}

template <typename RET>
static void
checkRandomFloatFloat (
    const std::string_view &funcname,
    const SampleRandomFloatFloatBase<RET> &funcs,
    const Description::Sample2Arg<typename RET::FloatType> &sample,
    const RoundSet &roundModes, FailMode failmode)
{
  using FloatType = typename RET::FloatType;

  std::vector<RngType> gens (rngStates.size ());

  for (auto &rnd : roundModes)
    {
      /* Reseed the RNG generator with the same seed to check the same numbers
	 for each rounding mode.  */
      for (unsigned i = 0; i < rngStates.size (); i++)
	gens[i] = RngType (rngStates[i]);

#pragma omp declare reduction(                                                \
	ulpAccumulatorReduction : UlpAccumulator<                             \
		FloatType> : ulpAccumulatorReduction(omp_out, omp_in))        \
    initializer(omp_priv = UlpAccumulator<FloatType> ())

      auto start = ClockType::now ();

      std::uniform_real_distribution<FloatType> distX (sample.arg_x.start,
						       sample.arg_x.end);
      std::uniform_real_distribution<FloatType> distY (sample.arg_y.start,
						       sample.arg_y.end);

      UlpAccumulator<FloatType> ulpaccrange;

#pragma omp parallel firstprivate(distX, distY, failmode) shared(sample, rnd)
      {
	RoundSetup<FloatType> roundSetup (rnd.mode);

#pragma omp for reduction(ulpAccumulatorReduction : ulpaccrange)
	for (std::uint64_t i = 0; i < sample.count; i++)
	  {
	    auto ret = funcs (gens[getThreadNum ()], distX, distY, rnd.mode);
	    if (!ret->check ())
	      switch (failmode)
		{
		case FailMode::FIRST:
#pragma omp critical
		  std::cerr << *ret;
		  std::exit (EXIT_FAILURE);
		  break;
		case FailMode::ALL:
#pragma omp critical
		  std::cerr << *ret;
		  [[fallthrough]];
		default:
		  break;
		}
	    ulpaccrange[ret->ulp] += 1;
	  }
      }

      printAccumulator (rnd.name, sample, ulpaccrange);

      auto end = ClockType::now ();
      printlnTimestamp (
	  "Elapsed time {}",
	  std::chrono::duration_cast<std::chrono::duration<double> > (
	      end - start));
      printlnTimestamp ("");
    }
}

template <typename RET>
static void
checkRandomFloatLLI (
    const std::string_view &funcname,
    const SampleRandomFloatLLIBase<RET> &funcs,
    const Description::Sample2ArgLli<typename RET::FloatType> &sample,
    const RoundSet &roundModes, FailMode failmode)
{
  using FloatType = typename RET::FloatType;
  using Arg2Type = long long int;

  std::vector<RngType> gens (rngStates.size ());

  for (auto &rnd : roundModes)
    {
      /* Reseed the RNG generator with the same seed to check the same numbers
	 for each rounding mode.  */
      for (unsigned i = 0; i < rngStates.size (); i++)
	gens[i] = RngType (rngStates[i]);

#pragma omp declare reduction(                                                \
	ulpAccumulatorReduction : UlpAccumulator<                             \
		FloatType> : ulpAccumulatorReduction(omp_out, omp_in))        \
    initializer(omp_priv = UlpAccumulator<FloatType> ())

      auto start = ClockType::now ();

      std::uniform_real_distribution<FloatType> distX (sample.arg_x.start,
						       sample.arg_x.end);
      std::uniform_int_distribution<Arg2Type> distY (sample.arg_y.start,
						     sample.arg_y.end);

      UlpAccumulator<FloatType> ulpaccrange;

#pragma omp parallel firstprivate(distX, distY, failmode) shared(sample, rnd)
      {
	RoundSetup<FloatType> roundSetup (rnd.mode);

#pragma omp for reduction(ulpAccumulatorReduction : ulpaccrange)
	for (std::uint64_t i = 0; i < sample.count; i++)
	  {
	    auto ret = funcs (gens[getThreadNum ()], distX, distY, rnd.mode);
	    if (!ret->check ())
	      switch (failmode)
		{
		case FailMode::FIRST:
#pragma omp critical
		  std::cerr << *ret;
		  std::exit (EXIT_FAILURE);
		  break;
		case FailMode::ALL:
#pragma omp critical
		  std::cerr << *ret;
		  [[fallthrough]];
		default:
		  break;
		}
	    ulpaccrange[ret->ulp] += 1;
	  }
      }

      printAccumulator (rnd.name, sample, ulpaccrange);

      auto end = ClockType::now ();
      printlnTimestamp (
	  "Elapsed time {}",
	  std::chrono::duration_cast<std::chrono::duration<double> > (
	      end - start));
      printlnTimestamp ("");
    }
}

template <typename RET>
static void
checkFull (const std::string_view &funcname, const SampleFull<RET> &funcs,
	   const Description::FullRange &sample, const RoundSet &roundModes,
	   FailMode failmode)
{
  using FloatType = typename RET::FloatType;

  for (auto &rnd : roundModes)
    {
#pragma omp declare reduction(                                                \
	ulpAccumulatorReduction : UlpAccumulator<                             \
		FloatType> : ulpAccumulatorReduction(omp_out, omp_in))        \
    initializer(omp_priv = UlpAccumulator<FloatType> ())

      UlpAccumulator<FloatType> ulpaccrange;

#pragma omp parallel firstprivate(failmode) shared(funcs, rnd)
      {
	RoundSetup<FloatType> roundSetup (rnd.mode);

// Out of range inputs might take way less time than normal one, also use
// a large chunk size to minimize the overhead from dynamic scheduline.
#pragma omp for reduction(ulpAccumulatorReduction : ulpaccrange)              \
    schedule(dynamic)
	for (std::uint64_t i = sample.start; i < sample.end; i++)
	  {
	    auto ret = funcs (i, rnd.mode);
	    if (!ret->checkFull ())
	      switch (failmode)
		{
		case FailMode::FIRST:
#pragma omp critical
		  std::cerr << *ret;
		  std::exit (EXIT_FAILURE);
		  break;
		case FailMode::ALL:
#pragma omp critical
		  std::cerr << *ret;
		  [[fallthrough]];
		default:
		  break;
		}
	    ulpaccrange[ret->ulp] += 1;
	  }
      }

      printAccumulator (rnd.name, sample, ulpaccrange);
      printlnTimestamp ("");
    }
}

template <typename F>
static void
runFloat (const Description &desc, const RoundSet &roundModes,
	  FailMode failmode, const std::string &max_ulp_str)
{
  auto func = getFunctionFloat<F> (desc.FunctionName).value ();
  if (!func.first)
    error ("libc does not provide {}", desc.FunctionName);

  const auto max_ulp = floatrange::fromStr<F> (max_ulp_str);
  if (!max_ulp)
    error ("invalid floating point: {}", max_ulp_str);

  printlnTimestamp ("Checking function {}", desc.FunctionName);
  printlnTimestamp ("");

  auto start = ClockType::now ();

  for (auto &sample : desc.Samples)
    {
      if (auto *psample = std::get_if<Description::Sample1Arg<F> > (&sample))
	checkRandomFloat (
	    desc.FunctionName,
	    RandomFloat<F>{ func.first, func.second, max_ulp.value () },
	    *psample, roundModes, failmode);
      else if (auto *psample = std::get_if<Description::FullRange> (&sample))
	checkFull (desc.FunctionName,
		   FullFloat<F>{ func.first, func.second, max_ulp.value () },
		   *psample, roundModes, failmode);
      else
	error ("invalid sample type");
    }

  auto end = ClockType::now ();
  printlnTimestamp (
      "Total elapsed time {}",
      std::chrono::duration_cast<std::chrono::duration<double> > (end
								  - start));
}

template <typename F>
static void
runFloatpFloatp (const Description &desc, const RoundSet &roundModes,
		 FailMode failmode, const std::string &max_ulp_str)
{
  auto func = getFunctionFloatpFloatp<F> (desc.FunctionName).value ();
  if (!func.first)
    error ("libc does not provide {}", desc.FunctionName);

  const auto max_ulp = floatrange::fromStr<F> (max_ulp_str);
  if (!max_ulp)
    error ("invalid floating point: {}", max_ulp_str);

  printlnTimestamp ("Checking function {}", desc.FunctionName);
  printlnTimestamp ("");

  auto start = ClockType::now ();

  for (auto &sample : desc.Samples)
    {
      if (auto *psample = std::get_if<Description::Sample1Arg<F> > (&sample))
	checkRandomFloat (
	    desc.FunctionName,
	    RandomFloatpFloatp<F>{ func.first, func.second, max_ulp.value () },
	    *psample, roundModes, failmode);
      else if (auto *psample = std::get_if<Description::FullRange> (&sample))
	checkFull (
	    desc.FunctionName,
	    FullFloatpFloatp<F>{ func.first, func.second, max_ulp.value () },
	    *psample, roundModes, failmode);
      else
	error ("invalid sample type");
    }

  auto end = ClockType::now ();
  printlnTimestamp (
      "Total elapsed time {}",
      std::chrono::duration_cast<std::chrono::duration<double> > (end
								  - start));
}

template <typename F>
static void
runFloatFloat (const Description &desc, const RoundSet &roundModes,
	       FailMode failmode, const std::string &max_ulp_str)
{
  auto func = getFunctionFloatFloat<F> (desc.FunctionName).value ();
  if (!func.first)
    error ("libc does not provide {}", desc.FunctionName);

  const auto max_ulp = floatrange::fromStr<F> (max_ulp_str);
  if (!max_ulp)
    error ("invalid floating point: {}", max_ulp_str);

  printlnTimestamp ("Checking FunctionName {}", desc.FunctionName);
  printlnTimestamp ("");

  auto start = ClockType::now ();

  for (auto &sample : desc.Samples)
    {
      if (auto *psample = std::get_if<Description::Sample2Arg<F> > (&sample))
	checkRandomFloatFloat (
	    desc.FunctionName,
	    RandomFloatFloat<F>{ func.first, func.second, max_ulp.value () },
	    *psample, roundModes, failmode);
      else
	error ("invalid sample type");
    }

  auto end = ClockType::now ();
  printlnTimestamp (
      "Total elapsed time {}",
      std::chrono::duration_cast<std::chrono::duration<double> > (end
								  - start));
}

template <typename F>
static void
runFloatLLI (const Description &desc, const RoundSet &roundModes,
	     FailMode failmode, const std::string &max_ulp_str)
{
  auto func = getFunctionFloatLLI<F> (desc.FunctionName).value ();
  if (!func.first)
    error ("libc does not provide {}", desc.FunctionName);

  const auto max_ulp = floatrange::fromStr<F> (max_ulp_str);
  if (!max_ulp)
    error ("invalid floating point: {}", max_ulp_str);

  printlnTimestamp ("Checking FunctionName {}", desc.FunctionName);
  printlnTimestamp ("");

  auto start = ClockType::now ();

  for (auto &sample : desc.Samples)
    {
      if (auto *psample
	  = std::get_if<Description::Sample2ArgLli<F> > (&sample))
	checkRandomFloatLLI (
	    desc.FunctionName,
	    RandomFloatLLI<F>{ func.first, func.second, max_ulp.value () },
	    *psample, roundModes, failmode);
      else
	error ("invalid sample type");
    }

  auto end = ClockType::now ();
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
      .default_value (defaultRoundOption ());

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

  const RoundSet roundModes
      = roundFromOption (options.get<std::string> ("-r"));

  FailMode failmode = failModeFromOptions (options.get<std::string> ("-f"));

  Description desc;
  if (auto r = desc.parse (options.get<std::string> ("-d")); !r)
    error ("{}", r.error ());

  std::string max_ulp = options.get<std::string> ("-m");

  initRandomState ();

  auto functype = getFunctionType (desc.FunctionName);
  if (!functype)
    error ("invalid FunctionName: {}", desc.FunctionName);

  switch (functype.value ())
    {
    case refimpls::FunctionType::f32_f:
      runFloat<float> (desc, roundModes, failmode, max_ulp);
      break;
    case refimpls::FunctionType::f64_f:
      runFloat<double> (desc, roundModes, failmode, max_ulp);
      break;

    case refimpls::FunctionType::f32_f_f:
      runFloatFloat<float> (desc, roundModes, failmode, max_ulp);
      break;
    case refimpls::FunctionType::f64_f_f:
      runFloatFloat<double> (desc, roundModes, failmode, max_ulp);
      break;

    case refimpls::FunctionType::f32_f_lli:
      runFloatLLI<float> (desc, roundModes, failmode, max_ulp);
      break;
    case refimpls::FunctionType::f64_f_lli:
      runFloatLLI<double> (desc, roundModes, failmode, max_ulp);
      break;

    case refimpls::FunctionType::f32_f_fp_fp:
      runFloatpFloatp<float> (desc, roundModes, failmode, max_ulp);
      break;
    case refimpls::FunctionType::f64_f_fp_fp:
      runFloatpFloatp<double> (desc, roundModes, failmode, max_ulp);
      break;

    default:
      error ("function type \"{}\" not implemented", functype.value ());
    }
}
