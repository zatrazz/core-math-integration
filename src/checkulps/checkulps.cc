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

#include <cerrno>
#include <fenv.h>
#include <omp.h>

#include "description.h"
#include "floatranges.h"
#include "iohelper.h"
#include "refimpls.h"
#include "wyhash64.h"
#include "strhelper.h"

// This is the threshold used by glibc that triggers a failure.
static constexpr auto kMaxUlpStr = "0.0";

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
      /* NaN/Inf have no meaningful ulp; return NaN so that ulpdiff() yields
	 NaN, which the Result constructor maps to a zero ulp distance.  */
      ulp = std::numeric_limits<F>::quiet_NaN ();
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

// ULP histogram: counts how many samples produced each ULP-distance value.
// The overwhelmingly common values are small non-negative integers (0, 1, 2,
// ...), so those are counted in a flat array on the hot path (a direct index,
// no allocation or tree walk); any fractional or large value falls back to a
// map.  This keeps the per-sample accumulation, which runs once per input per
// rounding mode, allocation-free in the common case.
template <typename F> class UlpAccumulator
{
public:
  static constexpr int kSmall = 64;

  inline void
  add (F ulp)
  {
    if (ulp >= 0 && ulp < static_cast<F> (kSmall))
      {
	int i = static_cast<int> (ulp);
	if (static_cast<F> (i) == ulp)
	  {
	    small_[i] += 1;
	    return;
	  }
      }
    rest_[ulp] += 1;
  }

  void
  merge (const UlpAccumulator &other)
  {
    for (int i = 0; i < kSmall; i++)
      small_[i] += other.small_[i];
    for (const auto &[k, v] : other.rest_)
      rest_[k] += v;
  }

  std::uint64_t
  total () const
  {
    std::uint64_t t = 0;
    for (auto c : small_)
      t += c;
    for (const auto &[k, v] : rest_)
      t += v;
    return t;
  }

  // Return the histogram as a sorted key -> count map, merging the small-value
  // array and the overflow map.  Only used for reporting, off the hot path.
  std::map<F, std::uint64_t>
  sorted () const
  {
    std::map<F, std::uint64_t> m (rest_);
    for (int i = 0; i < kSmall; i++)
      if (small_[i])
	m[static_cast<F> (i)] += small_[i];
    return m;
  }

private:
  std::array<std::uint64_t, kSmall> small_{};
  std::map<F, std::uint64_t> rest_;
};

// One ULP histogram per rounding mode.  The reference (MPFR) result is computed
// once per input and rounded to all selected modes, so the driver accumulates
// into a per-mode set indexed by REF_RND*.
template <typename F>
using UlpAccumulatorSet = std::array<UlpAccumulator<F>, REF_NRND>;

template <typename F>
static void
ulpAccumulatorSetReduction (UlpAccumulatorSet<F> &inout,
			    UlpAccumulatorSet<F> &in)
{
  for (int i = 0; i < REF_NRND; i++)
    inout[i].merge (in[i]);
}

// Map a C99 rounding mode (FE_*) to its reference result index (REF_RND*).
static int
refIndex (int mode)
{
  switch (mode)
    {
    case FE_TONEAREST:
      return REF_RNDN;
    case FE_UPWARD:
      return REF_RNDU;
    case FE_DOWNWARD:
      return REF_RNDD;
    case FE_TOWARDZERO:
      return REF_RNDZ;
    default:
      std::unreachable ();
    }
}

// Build the reference rounding-mode bitmask (1u << REF_RND*) selecting which
// modes the reference functions should round to.
static unsigned
maskFromRoundSet (const RoundSet &roundModes)
{
  unsigned mask = 0;
  for (const auto &rnd : roundModes)
    mask |= 1u << refIndex (rnd.mode);
  return mask;
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
  const std::uint64_t ulptotal = ulpacc.total ();

  printlnTimestamp (
      "Checking rounding mode {:13}, range [{:9.2g},{:9.2g}], count {}",
      rndname, sample.arg.start, sample.arg.end, ulptotal);

  for (const auto &ulp : ulpacc.sorted ())
    printlnTimestamp ("    {:g}: {:16} {:6.2f}%", ulp.first, ulp.second,
		      ((double) ulp.second / (double) ulptotal) * 100.0);
}

template <typename F>
static void
printAccumulator (const std::string_view &rndname,
		  const Description::Sample2Arg<F> &sample,
		  const UlpAccumulator<F> &ulpacc)
{
  const std::uint64_t ulptotal = ulpacc.total ();

  printlnTimestamp ("Checking rounding mode {:13}, range x=[{:9.2g},{:9.2g}], "
		    "y=[{:9.2g},{:9.2g}], count {}",
		    rndname, sample.arg_x.start, sample.arg_x.end,
		    sample.arg_y.start, sample.arg_y.end, ulptotal);

  for (const auto &ulp : ulpacc.sorted ())
    printlnTimestamp ("    {:g}: {:16} {:6.2f}%", ulp.first, ulp.second,
		      ((double) ulp.second / (double) ulptotal) * 100.0);
}

template <typename F>
static void
printAccumulator (const std::string_view &rndname,
		  const Description::Sample2ArgLli<F> &sample,
		  const UlpAccumulator<F> &ulpacc)
{
  const std::uint64_t ulptotal = ulpacc.total ();

  printlnTimestamp ("Checking rounding mode {:13}, range x=[{:9.2g},{:9.2g}], "
		    "y=[{},{}], count {}",
		    rndname, sample.arg_x.start, sample.arg_x.end,
		    sample.arg_y.start, sample.arg_y.end, ulptotal);

  for (const auto &ulp : ulpacc.sorted ())
    printlnTimestamp ("    {:g}: {:16} {:6.2f}%", ulp.first, ulp.second,
		      ((double) ulp.second / (double) ulptotal) * 100.0);
}

template <typename F>
static void
printAccumulator (const std::string_view &rndname,
		  const Description::FullRange &sample,
		  const UlpAccumulator<F> &ulpacc)
{
  const std::uint64_t ulptotal = ulpacc.total ();

  printlnTimestamp ("Checking rounding mode {:13}, {}", rndname, sample.name);

  for (const auto &ulp : ulpacc.sorted ())
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
  }

  bool
  check (void) const
  {
    return ulp <= max;
  }

  bool
  checkFull (void) const
  {
    if (isSignaling (computed) || isSignaling (expected))
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

  virtual void printTo (std::ostream &) const = 0;

  const RoundMode &roundMode;
  FloatType computed;
  FloatType expected;
  FloatType ulp;
  FloatType max;
};

template <typename F> struct std::formatter<Result<F> >
{
  constexpr auto
  parse (std::format_parse_context &ctx)
  {
    return ctx.begin ();
  }

  auto
  format (const Result<F> &obj, std::format_context &ctx) const
  {
    std::ostringstream oss;
    obj.printTo (oss);
    return std::format_to (ctx.out (), "{}", oss.str ());
  }
};

template <typename F> struct ResultFloat : public Result<F>
{
  explicit ResultFloat (int r, F i, F c, F e, F m)
      : Result<F> (r, c, e, m), input (i)
  {
  }

  void
  printTo (std::ostream &os) const override
  {
    os << std::format (
	"{} ulp={:1.0f} input={:#a} computed={:#a} expected={:#a}",
	Result<F>::roundMode.name, Result<F>::ulp, input, Result<F>::computed,
	Result<F>::expected);
  }

  F input;
};

template <typename F>
struct std::formatter<ResultFloat<F> > : std::formatter<Result<F> >
{
  // No body needed; it inherits 'parse' and 'format' from the base.
  // When formatting ResultFloat<F>, it is implicitly cast to Result<F>&
  // to match the signature of the base 'format' function.
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
    return ulp <= max;
  }

  bool
  checkFull (void) const
  {
    if (isSignaling (computed1) || isSignaling (computed2)
	|| isSignaling (expected2) || isSignaling (expected2))
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

  void
  printTo (std::ostream &os) const
  {
    os << std::format ("{} ulp={:1.0f} input={:#a} computed=({:#a} {:#a}) "
		       "expected=({:#a} {:#a})",
		       roundMode.name, ulp, input, computed1, computed2,
		       expected1, expected2);
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

template <typename F> struct std::formatter<ResultFloatpFloatp<F> >
{
  constexpr auto
  parse (std::format_parse_context &ctx)
  {
    return ctx.begin ();
  }

  auto
  format (const ResultFloatpFloatp<F> &obj, std::format_context &ctx) const
  {
    std::ostringstream oss;
    obj.printTo (oss);
    return std::format_to (ctx.out (), "{}", oss.str ());
  }
};

template <typename F> struct ResultFloatFloat : public Result<F>
{
  explicit ResultFloatFloat (int r, F i0, F i1, F c, F e, F m)
      : Result<F> (r, c, e, m), input0 (i0), input1 (i1)
  {
  }

  void
  printTo (std::ostream &os) const override
  {
    os << std::format ("{} ulp={:1.0f} input=({:#a},{:#a}) computed={:#a} "
		       "expected={:#a}",
		       Result<F>::roundMode.name, Result<F>::ulp, input0,
		       input1, Result<F>::computed, Result<F>::expected);
  }

  F input0;
  F input1;
};

template <typename F>
struct std::formatter<ResultFloatFloat<F> > : std::formatter<Result<F> >
{
  // No body needed; it inherits 'parse' and 'format' from the base.
  // When formatting ResultFloatFloat<F>, it is implicitly cast to Result<F>&
  // to match the signature of the base 'format' function.
};

template <typename F> struct ResultFloatLLI : public Result<F>
{
public:
  explicit ResultFloatLLI (int r, F i0, long long int i1, F c, F e, F m)
      : Result<F> (r, c, e, m), input0 (i0), input1 (i1)
  {
  }

  void
  printTo (std::ostream &os) const override
  {
    os << std::format ("{} ulp={:1.0f} input=({:#a},{}) computed={:#a} "
		       "expected={:#a}",
		       Result<F>::roundMode.name, Result<F>::ulp, input0,
		       input1, Result<F>::computed, Result<F>::expected);
  }

  F input0;
  long long int input1;
};

template <typename F>
struct std::formatter<ResultFloatLLI<F> > : std::formatter<Result<F> >
{
  // No body needed; it inherits 'parse' and 'format' from the base.
  // When formatting ResultFloatLLI<F>, it is implicitly cast to Result<F>&
  // to match the signature of the base 'format' function.
};

// Number of inputs processed per tile.  Within a tile the reference results
// are computed once for all inputs, then each rounding mode is swept in one
// pass so the hardware rounding mode (fesetround, which also stalls the FP
// pipeline) is set once per mode per tile instead of once per input.
static constexpr std::uint64_t kTileSize = 256;

// The floating-point exceptions compared between the tested function and the
// reference.  Set by the --exceptions option.
// The exceptions compared between the tested function and the reference.
// FE_INEXACT is deliberately excluded: like glibc's libm tests, which mark
// inexact optional for every non-exact function, the inexact flag is not
// checked (every function handled here is non-exact / transcendental).
static constexpr int kDriverExcMask
    = FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW;
static bool gCheckExc = false;
static bool gCheckErrno = false;
// True when either option needs the reference exception side-channel (errno
// expectations are derived from the expected exceptions).
static bool gComputeExc = false;

// The errno value a correctly-rounded result with exceptions EXC and value
// EXPECTEDVALUE is expected to set, following glibc: EDOM for a domain error
// (invalid) and ERANGE for a pole (divide-by-zero), for overflow only when the
// result is infinite, and for underflow only when the result flushes to zero.
// For overflow to the maximum finite value or a nonzero subnormal underflow,
// glibc leaves errno unchanged (it treats ERANGE as optional there).
template <typename F>
static inline int
expectedErrno (unsigned exc, F expectedValue)
{
  if (exc & FE_INVALID)
    return EDOM;
  if (exc & FE_DIVBYZERO)
    return ERANGE;
  if ((exc & FE_OVERFLOW) && std::isinf (expectedValue))
    return ERANGE;
  if ((exc & FE_UNDERFLOW) && expectedValue == static_cast<F> (0))
    return ERANGE;
  return 0;
}

static std::string
errnoToStr (int e)
{
  switch (e)
    {
    case 0:
      return "0";
    case EDOM:
      return "EDOM";
    case ERANGE:
      return "ERANGE";
    default:
      return std::to_string (e);
    }
}

// ULP distance with the same NaN/Inf handling as the Result constructors: a
// NaN/Inf ulpdiff (from a NaN/Inf operand) maps to a zero distance.
template <typename F>
static inline F
ulpDistance (F computed, F expected)
{
  F u = ulpdiff (computed, expected);
  if (std::isnan (u) || std::isinf (u))
    return 0.0;
  return u;
}

// As ulpDistance, but saturated at MAX_ULP, matching ResultFloatpFloatp.
template <typename F>
static inline F
ulpDistanceClamped (F computed, F expected, F max_ulp)
{
  F u = ulpDistance (computed, expected);
  return u >= max_ulp ? max_ulp : u;
}

// Render an FE_* exception mask as a human-readable string.
static std::string
excToStr (unsigned e)
{
  if ((e & kDriverExcMask) == 0)
    return "none";
  std::string s;
  auto add = [&] (int bit, const char *name) {
    if (e & bit)
      {
	if (!s.empty ())
	  s += '|';
	s += name;
      }
  };
  add (FE_INVALID, "invalid");
  add (FE_DIVBYZERO, "divbyzero");
  add (FE_OVERFLOW, "overflow");
  add (FE_UNDERFLOW, "underflow");
  add (FE_INEXACT, "inexact");
  return s;
}

// Copy the reference exception side-channel (published by the last reference
// evaluation on this thread) into DST, masked to the checked exceptions.
static inline void
captureExpExc (unsigned dst[REF_NRND])
{
  for (int i = 0; i < REF_NRND; i++)
    dst[i] = refimpls_last_exc[i] & kDriverExcMask;
}

// Print a failing result (and exit for FailMode::FIRST).  The Result object is
// only constructed on the (rare) failure path.
template <typename RET>
static inline void
reportFailure (const RET &ret, FailMode failmode)
{
#pragma omp critical
  {
    printlnErrorTimestamp ("{}", ret);
    if (failmode == FailMode::FIRST)
      std::exit (EXIT_FAILURE);
  }
}

// Print an exception mismatch (and exit for FailMode::FIRST).
template <typename RET>
static inline void
reportExcMismatch (const RET &ret, unsigned expected, unsigned raised,
		   FailMode failmode)
{
#pragma omp critical
  {
    printlnErrorTimestamp ("{} fexcept expected={} raised={}", ret,
			   excToStr (expected), excToStr (raised));
    if (failmode == FailMode::FIRST)
      std::exit (EXIT_FAILURE);
  }
}

// Print an errno mismatch (and exit for FailMode::FIRST).
template <typename RET>
static inline void
reportErrnoMismatch (const RET &ret, int expected, int got, FailMode failmode)
{
#pragma omp critical
  {
    printlnErrorTimestamp ("{} errno expected={} got={}", ret,
			   errnoToStr (expected), errnoToStr (got));
    if (failmode == FailMode::FIRST)
      std::exit (EXIT_FAILURE);
  }
}

// Record one sample into the histogram, reporting a ULP failure (when
// REPORT_VALUE), an exception mismatch (when DO_EXC) and/or an errno mismatch
// (when DO_ERRNO).  MK is a factory that builds the printable Result, invoked
// only on the failure path.
template <typename F, typename MakeRet>
static inline void
recordSample (F u, F max_ulp, bool reportValue, unsigned raised,
	      unsigned expExc, F expVal, bool doExc, int gotErrno,
	      bool doErrno, FailMode failmode, UlpAccumulator<F> &acc,
	      MakeRet mk)
{
  if (reportValue && u > max_ulp && failmode != FailMode::NONE)
    reportFailure (mk (), failmode);
  if (doExc && raised != expExc)
    reportExcMismatch (mk (), expExc, raised, failmode);
  if (doErrno)
    {
      int expErr = expectedErrno (expExc, expVal);
      if (gotErrno != expErr)
	reportErrnoMismatch (mk (), expErr, gotErrno, failmode);
    }
  acc.add (u);
}

//
// Random-sampling checks.  For each input the reference (MPFR) result is
// computed once and rounded to every selected mode; inputs are processed in
// tiles so the tested libc function is evaluated for a whole tile under each
// hardware rounding mode before switching to the next.
//

template <typename F>
static void
checkRandomFloat (const std::string_view &funcname, FuncF<F> func,
		  const FuncFReference<F> &ref, F max_ulp,
		  const Description::Sample1Arg<F> &sample,
		  const RoundSet &roundModes, FailMode failmode)
{
  using FloatType = F;
  const unsigned mask = maskFromRoundSet (roundModes);

  refimpls::setupReferenceImpl<FloatType> ();

  std::vector<RngType> gens (rngStates.size ());
  for (unsigned i = 0; i < rngStates.size (); i++)
    gens[i] = RngType (rngStates[i]);

  auto start = ClockType::now ();

  std::uniform_real_distribution<FloatType> dist (sample.arg.start,
						  sample.arg.end);

  UlpAccumulatorSet<FloatType> ulpacc;
  const std::uint64_t count = sample.count;
  const std::uint64_t ntiles = (count + kTileSize - 1) / kTileSize;

#pragma omp declare reduction(                                                \
      ulpAccumulatorSetReduction : UlpAccumulatorSet<                         \
	  FloatType> : ulpAccumulatorSetReduction(omp_out, omp_in))           \
      initializer(omp_priv = UlpAccumulatorSet<FloatType> ())

#pragma omp parallel firstprivate(dist, failmode) shared(roundModes)
  {
    int savedRound = fegetround ();
    std::array<FloatType, kTileSize> inbuf;
    std::array<std::array<FloatType, REF_NRND>, kTileSize> expbuf;
    std::array<std::array<unsigned, REF_NRND>, kTileSize> expexc;

#pragma omp for reduction(ulpAccumulatorSetReduction : ulpacc)
    for (std::uint64_t t = 0; t < ntiles; t++)
      {
	std::uint64_t base = t * kTileSize;
	std::uint64_t n = std::min<std::uint64_t> (kTileSize, count - base);
	RngType &gen = gens[getThreadNum ()];

	for (std::uint64_t j = 0; j < n; j++)
	  {
	    inbuf[j] = dist (gen);
	    ref (inbuf[j], mask, expbuf[j].data ());
	    if (gComputeExc)
	      captureExpExc (expexc[j].data ());
	  }

	for (const auto &rnd : roundModes)
	  {
	    fesetround (rnd.mode);
	    int idx = refIndex (rnd.mode);
	    UlpAccumulator<FloatType> &acc = ulpacc[idx];
	    for (std::uint64_t j = 0; j < n; j++)
	      {
		if (gCheckExc)
		  feclearexcept (kDriverExcMask);
		if (gCheckErrno)
		  errno = 0;
		FloatType computed = func (inbuf[j]);
		unsigned raised
		    = gCheckExc ? (unsigned) fetestexcept (kDriverExcMask) : 0u;
		int gotErrno = gCheckErrno ? errno : 0;
		FloatType u = ulpDistance (computed, expbuf[j][idx]);
		recordSample (u, max_ulp, true, raised, expexc[j][idx], expbuf[j][idx],
			      gCheckExc, gotErrno,
			      gCheckErrno, failmode, acc, [&] {
				return ResultFloat<FloatType> (
				    rnd.mode, inbuf[j], computed,
				    expbuf[j][idx], max_ulp);
			      });
	      }
	  }
      }

    fesetround (savedRound);
  }

  for (const auto &rnd : roundModes)
    printAccumulator (rnd.name, sample, ulpacc[refIndex (rnd.mode)]);

  auto end = ClockType::now ();
  printlnTimestamp (
      "Elapsed time {}",
      std::chrono::duration_cast<std::chrono::duration<double> > (end
								  - start));
  printlnTimestamp ("");
}

template <typename F>
static void
checkRandomFloatpFloatp (const std::string_view &funcname, FuncFpFp<F> func,
			 const FuncFpFpReference<F> &ref, F max_ulp,
			 const Description::Sample1Arg<F> &sample,
			 const RoundSet &roundModes, FailMode failmode)
{
  using FloatType = F;
  const unsigned mask = maskFromRoundSet (roundModes);

  refimpls::setupReferenceImpl<FloatType> ();

  std::vector<RngType> gens (rngStates.size ());
  for (unsigned i = 0; i < rngStates.size (); i++)
    gens[i] = RngType (rngStates[i]);

  auto start = ClockType::now ();

  std::uniform_real_distribution<FloatType> dist (sample.arg.start,
						  sample.arg.end);

  UlpAccumulatorSet<FloatType> ulpacc;
  const std::uint64_t count = sample.count;
  const std::uint64_t ntiles = (count + kTileSize - 1) / kTileSize;

#pragma omp declare reduction(                                                \
      ulpAccumulatorSetReduction : UlpAccumulatorSet<                         \
	  FloatType> : ulpAccumulatorSetReduction(omp_out, omp_in))           \
      initializer(omp_priv = UlpAccumulatorSet<FloatType> ())

#pragma omp parallel firstprivate(dist, failmode) shared(roundModes)
  {
    int savedRound = fegetround ();
    std::array<FloatType, kTileSize> inbuf;
    std::array<std::array<FloatType, REF_NRND>, kTileSize> expbuf0, expbuf1;
    std::array<std::array<unsigned, REF_NRND>, kTileSize> expexc;

#pragma omp for reduction(ulpAccumulatorSetReduction : ulpacc)
    for (std::uint64_t t = 0; t < ntiles; t++)
      {
	std::uint64_t base = t * kTileSize;
	std::uint64_t n = std::min<std::uint64_t> (kTileSize, count - base);
	RngType &gen = gens[getThreadNum ()];

	for (std::uint64_t j = 0; j < n; j++)
	  {
	    inbuf[j] = dist (gen);
	    ref (inbuf[j], mask, expbuf0[j].data (), expbuf1[j].data ());
	    if (gComputeExc)
	      captureExpExc (expexc[j].data ());
	  }

	for (const auto &rnd : roundModes)
	  {
	    fesetround (rnd.mode);
	    int idx = refIndex (rnd.mode);
	    UlpAccumulator<FloatType> &acc = ulpacc[idx];
	    for (std::uint64_t j = 0; j < n; j++)
	      {
		if (gCheckExc)
		  feclearexcept (kDriverExcMask);
		if (gCheckErrno)
		  errno = 0;
		FloatType computed0, computed1;
		func (inbuf[j], &computed0, &computed1);
		unsigned raised
		    = gCheckExc ? (unsigned) fetestexcept (kDriverExcMask) : 0u;
		int gotErrno = gCheckErrno ? errno : 0;
		FloatType u0
		    = ulpDistanceClamped (computed0, expbuf0[j][idx], max_ulp);
		FloatType u1
		    = ulpDistanceClamped (computed1, expbuf1[j][idx], max_ulp);
		FloatType u = u1 > u0 ? u1 : u0;
		recordSample (u, max_ulp, true, raised, expexc[j][idx], expbuf0[j][idx],
			      gCheckExc, gotErrno,
			      gCheckErrno, failmode, acc, [&] {
				return ResultFloatpFloatp<FloatType> (
				    rnd.mode, inbuf[j], computed0, computed1,
				    expbuf0[j][idx], expbuf1[j][idx], max_ulp);
			      });
	      }
	  }
      }

    fesetround (savedRound);
  }

  for (const auto &rnd : roundModes)
    printAccumulator (rnd.name, sample, ulpacc[refIndex (rnd.mode)]);

  auto end = ClockType::now ();
  printlnTimestamp (
      "Elapsed time {}",
      std::chrono::duration_cast<std::chrono::duration<double> > (end
								  - start));
  printlnTimestamp ("");
}

template <typename F>
static void
checkRandomFloatFloat (const std::string_view &funcname, FuncFF<F> func,
		       const FuncFFReference<F> &ref, F max_ulp,
		       const Description::Sample2Arg<F> &sample,
		       const RoundSet &roundModes, FailMode failmode)
{
  using FloatType = F;
  const unsigned mask = maskFromRoundSet (roundModes);

  refimpls::setupReferenceImpl<FloatType> ();

  std::vector<RngType> gens (rngStates.size ());
  for (unsigned i = 0; i < rngStates.size (); i++)
    gens[i] = RngType (rngStates[i]);

  auto start = ClockType::now ();

  std::uniform_real_distribution<FloatType> distX (sample.arg_x.start,
						   sample.arg_x.end);
  std::uniform_real_distribution<FloatType> distY (sample.arg_y.start,
						   sample.arg_y.end);

  UlpAccumulatorSet<FloatType> ulpacc;
  const std::uint64_t count = sample.count;
  const std::uint64_t ntiles = (count + kTileSize - 1) / kTileSize;

#pragma omp declare reduction(                                                \
      ulpAccumulatorSetReduction : UlpAccumulatorSet<                         \
	  FloatType> : ulpAccumulatorSetReduction(omp_out, omp_in))           \
      initializer(omp_priv = UlpAccumulatorSet<FloatType> ())

#pragma omp parallel firstprivate(distX, distY, failmode) shared(roundModes)
  {
    int savedRound = fegetround ();
    std::array<FloatType, kTileSize> inbuf0, inbuf1;
    std::array<std::array<FloatType, REF_NRND>, kTileSize> expbuf;
    std::array<std::array<unsigned, REF_NRND>, kTileSize> expexc;

#pragma omp for reduction(ulpAccumulatorSetReduction : ulpacc)
    for (std::uint64_t t = 0; t < ntiles; t++)
      {
	std::uint64_t base = t * kTileSize;
	std::uint64_t n = std::min<std::uint64_t> (kTileSize, count - base);
	RngType &gen = gens[getThreadNum ()];

	for (std::uint64_t j = 0; j < n; j++)
	  {
	    inbuf0[j] = distX (gen);
	    inbuf1[j] = distY (gen);
	    ref (inbuf0[j], inbuf1[j], mask, expbuf[j].data ());
	    if (gComputeExc)
	      captureExpExc (expexc[j].data ());
	  }

	for (const auto &rnd : roundModes)
	  {
	    fesetround (rnd.mode);
	    int idx = refIndex (rnd.mode);
	    UlpAccumulator<FloatType> &acc = ulpacc[idx];
	    for (std::uint64_t j = 0; j < n; j++)
	      {
		if (gCheckExc)
		  feclearexcept (kDriverExcMask);
		if (gCheckErrno)
		  errno = 0;
		FloatType computed = func (inbuf0[j], inbuf1[j]);
		unsigned raised
		    = gCheckExc ? (unsigned) fetestexcept (kDriverExcMask) : 0u;
		int gotErrno = gCheckErrno ? errno : 0;
		FloatType u = ulpDistance (computed, expbuf[j][idx]);
		recordSample (u, max_ulp, true, raised, expexc[j][idx], expbuf[j][idx],
			      gCheckExc, gotErrno,
			      gCheckErrno, failmode, acc, [&] {
				return ResultFloatFloat<FloatType> (
				    rnd.mode, inbuf0[j], inbuf1[j], computed,
				    expbuf[j][idx], max_ulp);
			      });
	      }
	  }
      }

    fesetround (savedRound);
  }

  for (const auto &rnd : roundModes)
    printAccumulator (rnd.name, sample, ulpacc[refIndex (rnd.mode)]);

  auto end = ClockType::now ();
  printlnTimestamp (
      "Elapsed time {}",
      std::chrono::duration_cast<std::chrono::duration<double> > (end
								  - start));
  printlnTimestamp ("");
}

template <typename F>
static void
checkRandomFloatLLI (const std::string_view &funcname, FuncFLLI<F> func,
		     const FuncFLLIReference<F> &ref, F max_ulp,
		     const Description::Sample2ArgLli<F> &sample,
		     const RoundSet &roundModes, FailMode failmode)
{
  using FloatType = F;
  using Arg2Type = long long int;
  const unsigned mask = maskFromRoundSet (roundModes);

  refimpls::setupReferenceImpl<FloatType> ();

  std::vector<RngType> gens (rngStates.size ());
  for (unsigned i = 0; i < rngStates.size (); i++)
    gens[i] = RngType (rngStates[i]);

  auto start = ClockType::now ();

  std::uniform_real_distribution<FloatType> distX (sample.arg_x.start,
						   sample.arg_x.end);
  std::uniform_int_distribution<Arg2Type> distY (sample.arg_y.start,
						 sample.arg_y.end);

  UlpAccumulatorSet<FloatType> ulpacc;
  const std::uint64_t count = sample.count;
  const std::uint64_t ntiles = (count + kTileSize - 1) / kTileSize;

#pragma omp declare reduction(                                                \
      ulpAccumulatorSetReduction : UlpAccumulatorSet<                         \
	  FloatType> : ulpAccumulatorSetReduction(omp_out, omp_in))           \
      initializer(omp_priv = UlpAccumulatorSet<FloatType> ())

#pragma omp parallel firstprivate(distX, distY, failmode) shared(roundModes)
  {
    int savedRound = fegetround ();
    std::array<FloatType, kTileSize> inbuf0;
    std::array<Arg2Type, kTileSize> inbuf1;
    std::array<std::array<FloatType, REF_NRND>, kTileSize> expbuf;
    std::array<std::array<unsigned, REF_NRND>, kTileSize> expexc;

#pragma omp for reduction(ulpAccumulatorSetReduction : ulpacc)
    for (std::uint64_t t = 0; t < ntiles; t++)
      {
	std::uint64_t base = t * kTileSize;
	std::uint64_t n = std::min<std::uint64_t> (kTileSize, count - base);
	RngType &gen = gens[getThreadNum ()];

	for (std::uint64_t j = 0; j < n; j++)
	  {
	    inbuf0[j] = distX (gen);
	    inbuf1[j] = distY (gen);
	    ref (inbuf0[j], inbuf1[j], mask, expbuf[j].data ());
	    if (gComputeExc)
	      captureExpExc (expexc[j].data ());
	  }

	for (const auto &rnd : roundModes)
	  {
	    fesetround (rnd.mode);
	    int idx = refIndex (rnd.mode);
	    UlpAccumulator<FloatType> &acc = ulpacc[idx];
	    for (std::uint64_t j = 0; j < n; j++)
	      {
		if (gCheckExc)
		  feclearexcept (kDriverExcMask);
		if (gCheckErrno)
		  errno = 0;
		FloatType computed = func (inbuf0[j], inbuf1[j]);
		unsigned raised
		    = gCheckExc ? (unsigned) fetestexcept (kDriverExcMask) : 0u;
		int gotErrno = gCheckErrno ? errno : 0;
		FloatType u = ulpDistance (computed, expbuf[j][idx]);
		recordSample (u, max_ulp, true, raised, expexc[j][idx], expbuf[j][idx],
			      gCheckExc, gotErrno,
			      gCheckErrno, failmode, acc, [&] {
				return ResultFloatLLI<FloatType> (
				    rnd.mode, inbuf0[j], inbuf1[j], computed,
				    expbuf[j][idx], max_ulp);
			      });
	      }
	  }
      }

    fesetround (savedRound);
  }

  for (const auto &rnd : roundModes)
    printAccumulator (rnd.name, sample, ulpacc[refIndex (rnd.mode)]);

  auto end = ClockType::now ();
  printlnTimestamp (
      "Elapsed time {}",
      std::chrono::duration_cast<std::chrono::duration<double> > (end
								  - start));
  printlnTimestamp ("");
}

//
// Full-range checks: iterate every input bit pattern in [start, end).  ULP
// failure reporting is intentionally disabled here, as the output would be
// dominated by out-of-domain inputs; only the ULP distribution is accumulated
// (and exception mismatches, when enabled, are reported).
//

template <typename F>
static void
checkFull (const std::string_view &funcname, FuncF<F> func,
	   const FuncFReference<F> &ref, F max_ulp,
	   const Description::FullRange &sample, const RoundSet &roundModes,
	   FailMode failmode)
{
  using FloatType = F;
  const unsigned mask = maskFromRoundSet (roundModes);

  refimpls::setupReferenceImpl<FloatType> ();

  UlpAccumulatorSet<FloatType> ulpacc;
  const std::uint64_t ntiles
      = (sample.end - sample.start + kTileSize - 1) / kTileSize;

#pragma omp declare reduction(                                                \
      ulpAccumulatorSetReduction : UlpAccumulatorSet<                         \
	  FloatType> : ulpAccumulatorSetReduction(omp_out, omp_in))           \
      initializer(omp_priv = UlpAccumulatorSet<FloatType> ())

#pragma omp parallel firstprivate(failmode) shared(roundModes)
  {
    int savedRound = fegetround ();
    std::array<FloatType, kTileSize> inbuf;
    std::array<std::array<FloatType, REF_NRND>, kTileSize> expbuf;
    std::array<std::array<unsigned, REF_NRND>, kTileSize> expexc;

// Out of range inputs might take way less time than normal ones; use dynamic
// scheduling so tiles are balanced across threads.
#pragma omp for reduction(ulpAccumulatorSetReduction : ulpacc) schedule(dynamic)
    for (std::uint64_t t = 0; t < ntiles; t++)
      {
	std::uint64_t base = sample.start + t * kTileSize;
	std::uint64_t n
	    = std::min<std::uint64_t> (kTileSize, sample.end - base);

	for (std::uint64_t j = 0; j < n; j++)
	  {
	    inbuf[j] = floatrange::Limits<FloatType>::from (base + j);
	    ref (inbuf[j], mask, expbuf[j].data ());
	    if (gComputeExc)
	      captureExpExc (expexc[j].data ());
	  }

	for (const auto &rnd : roundModes)
	  {
	    fesetround (rnd.mode);
	    int idx = refIndex (rnd.mode);
	    UlpAccumulator<FloatType> &acc = ulpacc[idx];
	    for (std::uint64_t j = 0; j < n; j++)
	      {
		if (gCheckExc)
		  feclearexcept (kDriverExcMask);
		if (gCheckErrno)
		  errno = 0;
		FloatType computed = func (inbuf[j]);
		unsigned raised
		    = gCheckExc ? (unsigned) fetestexcept (kDriverExcMask) : 0u;
		int gotErrno = gCheckErrno ? errno : 0;
		FloatType u = ulpDistance (computed, expbuf[j][idx]);
		recordSample (u, max_ulp, false, raised, expexc[j][idx], expbuf[j][idx],
			      gCheckExc, gotErrno,
			      gCheckErrno, failmode, acc, [&] {
				return ResultFloat<FloatType> (
				    rnd.mode, inbuf[j], computed,
				    expbuf[j][idx], max_ulp);
			      });
	      }
	  }
      }

    fesetround (savedRound);
  }

  for (const auto &rnd : roundModes)
    printAccumulator (rnd.name, sample, ulpacc[refIndex (rnd.mode)]);
  printlnTimestamp ("");
}

template <typename F>
static void
checkFullFloatpFloatp (const std::string_view &funcname, FuncFpFp<F> func,
		       const FuncFpFpReference<F> &ref, F max_ulp,
		       const Description::FullRange &sample,
		       const RoundSet &roundModes, FailMode failmode)
{
  using FloatType = F;
  const unsigned mask = maskFromRoundSet (roundModes);

  refimpls::setupReferenceImpl<FloatType> ();

  UlpAccumulatorSet<FloatType> ulpacc;
  const std::uint64_t ntiles
      = (sample.end - sample.start + kTileSize - 1) / kTileSize;

#pragma omp declare reduction(                                                \
      ulpAccumulatorSetReduction : UlpAccumulatorSet<                         \
	  FloatType> : ulpAccumulatorSetReduction(omp_out, omp_in))           \
      initializer(omp_priv = UlpAccumulatorSet<FloatType> ())

#pragma omp parallel firstprivate(failmode) shared(roundModes)
  {
    int savedRound = fegetround ();
    std::array<FloatType, kTileSize> inbuf;
    std::array<std::array<FloatType, REF_NRND>, kTileSize> expbuf0, expbuf1;
    std::array<std::array<unsigned, REF_NRND>, kTileSize> expexc;

#pragma omp for reduction(ulpAccumulatorSetReduction : ulpacc) schedule(dynamic)
    for (std::uint64_t t = 0; t < ntiles; t++)
      {
	std::uint64_t base = sample.start + t * kTileSize;
	std::uint64_t n
	    = std::min<std::uint64_t> (kTileSize, sample.end - base);

	for (std::uint64_t j = 0; j < n; j++)
	  {
	    inbuf[j] = floatrange::Limits<FloatType>::from (base + j);
	    ref (inbuf[j], mask, expbuf0[j].data (), expbuf1[j].data ());
	    if (gComputeExc)
	      captureExpExc (expexc[j].data ());
	  }

	for (const auto &rnd : roundModes)
	  {
	    fesetround (rnd.mode);
	    int idx = refIndex (rnd.mode);
	    UlpAccumulator<FloatType> &acc = ulpacc[idx];
	    for (std::uint64_t j = 0; j < n; j++)
	      {
		if (gCheckExc)
		  feclearexcept (kDriverExcMask);
		if (gCheckErrno)
		  errno = 0;
		FloatType computed0, computed1;
		func (inbuf[j], &computed0, &computed1);
		unsigned raised
		    = gCheckExc ? (unsigned) fetestexcept (kDriverExcMask) : 0u;
		int gotErrno = gCheckErrno ? errno : 0;
		FloatType u0
		    = ulpDistanceClamped (computed0, expbuf0[j][idx], max_ulp);
		FloatType u1
		    = ulpDistanceClamped (computed1, expbuf1[j][idx], max_ulp);
		FloatType u = u1 > u0 ? u1 : u0;
		recordSample (u, max_ulp, false, raised, expexc[j][idx], expbuf0[j][idx],
			      gCheckExc, gotErrno,
			      gCheckErrno, failmode, acc, [&] {
				return ResultFloatpFloatp<FloatType> (
				    rnd.mode, inbuf[j], computed0, computed1,
				    expbuf0[j][idx], expbuf1[j][idx], max_ulp);
			      });
	      }
	  }
      }

    fesetround (savedRound);
  }

  for (const auto &rnd : roundModes)
    printAccumulator (rnd.name, sample, ulpacc[refIndex (rnd.mode)]);
  printlnTimestamp ("");
}

//
// Explicit value-list check: prints every result grouped by rounding mode.
// The reference result for each value is computed once and cached.
//

template <typename F>
static void
checkList (const std::string_view &funcname, const std::vector<F> &values,
	   FuncF<F> func, const FuncFReference<F> &ref, F max_ulp,
	   const RoundSet &roundModes, FailMode failmode)
{
  const unsigned mask = maskFromRoundSet (roundModes);

  // Compute the reference result (and expected exceptions) for every value
  // once, up front.
  std::vector<std::array<F, REF_NRND> > expected (values.size ());
  std::vector<std::array<unsigned, REF_NRND> > expexc (values.size ());
  {
    RoundSetup<F> roundSetup (FE_TONEAREST);
    for (std::size_t i = 0; i < values.size (); i++)
      {
	ref (values[i], mask, expected[i].data ());
	if (gComputeExc)
	  captureExpExc (expexc[i].data ());
      }
  }

  for (const auto &rnd : roundModes)
    {
      int idx = refIndex (rnd.mode);
      RoundSetup<F> roundSetup (rnd.mode);

      for (std::size_t i = 0; i < values.size (); i++)
	{
	  if (gCheckExc)
	    feclearexcept (kDriverExcMask);
	  if (gCheckErrno)
	    errno = 0;
	  F computed = func (values[i]);
	  unsigned raised
	      = gCheckExc ? (unsigned) fetestexcept (kDriverExcMask) : 0u;
	  int gotErrno = gCheckErrno ? errno : 0;
	  ResultFloat<F> ret (rnd.mode, values[i], computed, expected[i][idx],
			      max_ulp);
	  if (!ret.checkFull ())
	    {
	      switch (failmode)
		{
		case FailMode::FIRST:
		case FailMode::ALL:
		  printlnErrorTimestamp ("{}", ret);
		  if (failmode == FailMode::FIRST)
		    std::exit (EXIT_FAILURE);
		  [[fallthrough]];
		default:
		  break;
		}
	    }
	  else
	    printlnTimestamp ("{}", ret);
	  if (gCheckExc && raised != expexc[i][idx])
	    reportExcMismatch (ret, expexc[i][idx], raised, failmode);
	  if (gCheckErrno)
	    {
	      int expErr = expectedErrno (expexc[i][idx], expected[i][idx]);
	      if (gotErrno != expErr)
		reportErrnoMismatch (ret, expErr, gotErrno, failmode);
	    }
	}
    }

  printlnTimestamp ("");
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
	checkRandomFloat (desc.FunctionName, func.first, func.second,
			  max_ulp.value (), *psample, roundModes, failmode);
      else if (auto *psample = std::get_if<Description::FullRange> (&sample))
	checkFull (desc.FunctionName, func.first, func.second, max_ulp.value (),
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
	checkRandomFloatpFloatp (desc.FunctionName, func.first, func.second,
				 max_ulp.value (), *psample, roundModes,
				 failmode);
      else if (auto *psample = std::get_if<Description::FullRange> (&sample))
	checkFullFloatpFloatp (desc.FunctionName, func.first, func.second,
			       max_ulp.value (), *psample, roundModes,
			       failmode);
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
	checkRandomFloatFloat (desc.FunctionName, func.first, func.second,
			       max_ulp.value (), *psample, roundModes,
			       failmode);
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
	checkRandomFloatLLI (desc.FunctionName, func.first, func.second,
			     max_ulp.value (), *psample, roundModes, failmode);
      else
	error ("invalid sample type");
    }

  auto end = ClockType::now ();
  printlnTimestamp (
      "Total elapsed time {}",
      std::chrono::duration_cast<std::chrono::duration<double> > (end
								  - start));
}

static void
handleDescription (const std::string &descFile, const RoundSet &roundModes,
		   FailMode failmode, const std::string &maxUlp)
{
  Description desc;
  if (auto r = desc.parse (descFile); !r)
    error ("{}", r.error ());

  initRandomState ();

  auto functype = getFunctionType (desc.FunctionName);
  if (!functype)
    error ("invalid FunctionName: {}", desc.FunctionName);

  switch (functype.value ())
    {
    case refimpls::FunctionType::f32_f:
      runFloat<float> (desc, roundModes, failmode, maxUlp);
      break;
    case refimpls::FunctionType::f64_f:
      runFloat<double> (desc, roundModes, failmode, maxUlp);
      break;

    case refimpls::FunctionType::f32_f_f:
      runFloatFloat<float> (desc, roundModes, failmode, maxUlp);
      break;
    case refimpls::FunctionType::f64_f_f:
      runFloatFloat<double> (desc, roundModes, failmode, maxUlp);
      break;

    case refimpls::FunctionType::f32_f_lli:
      runFloatLLI<float> (desc, roundModes, failmode, maxUlp);
      break;
    case refimpls::FunctionType::f64_f_lli:
      runFloatLLI<double> (desc, roundModes, failmode, maxUlp);
      break;

    case refimpls::FunctionType::f32_f_fp_fp:
      runFloatpFloatp<float> (desc, roundModes, failmode, maxUlp);
      break;
    case refimpls::FunctionType::f64_f_fp_fp:
      runFloatpFloatp<double> (desc, roundModes, failmode, maxUlp);
      break;

    default:
      error ("function type \"{}\" not implemented", functype.value ());
    }
}

template <typename F>
static void
runFloatList (const std::string &functionName, const std::vector<F> &values,
	      const RoundSet &roundModes, FailMode failmode,
	      const std::string &maxUlpStr)
{
  auto func = getFunctionFloat<F> (functionName).value ();
  if (!func.first)
    error ("libc does not provide {}", functionName);

  const auto maxUlp = floatrange::fromStr<F> (maxUlpStr);
  if (!maxUlp)
    error ("invalid floating point: {}", maxUlpStr);

  printlnTimestamp ("Checking function {}", functionName);
  printlnTimestamp ("");

  auto start = ClockType::now ();

  checkList (functionName, values, func.first, func.second, maxUlp.value (),
	     roundModes, failmode);

  auto end = ClockType::now ();
  printlnTimestamp (
      "Total elapsed time {}",
      std::chrono::duration_cast<std::chrono::duration<double> > (end
								  - start));
}

template <typename F>
static std::vector<F>
stringListToFPList (const std::vector<std::string> &valueList)
{
  auto view = valueList | std::views::transform ([] (const std::string &s) {
		auto result = floatrange::fromStr<F> (s);
		if (!result.has_value ())
		  error ("invalid number: {}", s);
		return result.value ();
	      });

  std::vector<F> numbers (view.begin (), view.end ());
  return numbers;
}

static void
handleList (const std::string &functionName,
	    const std::vector<std::string> &values, const RoundSet &roundModes,
	    FailMode failmode, const std::string &maxUlp)
{
  auto functype = getFunctionType (functionName);
  if (!functype)
    error ("invalid FunctionName: {}", functionName);

  switch (functype.value ())
    {
    case refimpls::FunctionType::f32_f:
      runFloatList<float> (functionName, stringListToFPList<float> (values),
			   roundModes, failmode, maxUlp);
      break;
    case refimpls::FunctionType::f64_f:
      runFloatList<double> (functionName, stringListToFPList<double> (values),
			    roundModes, failmode, maxUlp);
      break;

#if 0
    case refimpls::FunctionType::f32_f_f:
      runFloatFloatValue<float> (desc, roundModes, failmode, maxUlp);
      break;
    case refimpls::FunctionType::f64_f_f:
      runFloatFloatValue<double> (desc, roundModes, failmode, maxUlp);
      break;

    case refimpls::FunctionType::f32_f_lli:
      runFloatLLIValue<float> (desc, roundModes, failmode, maxUlp);
      break;
    case refimpls::FunctionType::f64_f_lli:
      runFloatLLIValue<double> (desc, roundModes, failmode, maxUlp);
      break;

    case refimpls::FunctionType::f32_f_fp_fp:
      runFloatpFloatpValue<float> (desc, roundModes, failmode, maxUlp);
      break;
    case refimpls::FunctionType::f64_f_fp_fp:
      runFloatpFloatpValue<double> (desc, roundModes, failmode, maxUlp);
      break;
#endif

    default:
      error ("function type \"{}\" not implemented", functype.value ());
    }
}

int
main (int argc, char *argv[])
{
  argparse::ArgumentParser options ("checkulps");

  options.add_argument ("--description", "-d")
      .help ("input JSON descriptiorn file");

  options.add_argument ("--symbol", "-s")
      .help ("math function to check")
      .nargs (1);

  options.add_argument ("--rounding", "-r")
      .help ("rounding modes to test")
      .default_value (defaultRoundOption ());

  options.add_argument ("--failure", "-f")
      .help ("failure mode")
      .default_value ("none");

  options.add_argument ("--maxulps", "-m")
      .help ("max ULP used in check")
      .default_value (kMaxUlpStr);

  options.add_argument ("--exceptions", "-e")
      .help ("also check the floating-point exceptions raised by the function "
	     "(invalid, divide-by-zero, overflow, underflow; inexact is not "
	     "checked, as in glibc's tests)")
      .flag ();

  options.add_argument ("--errno", "-E")
      .help ("also check errno (EDOM/ERANGE) set by the function")
      .flag ();

  options.add_argument ("values")
      .nargs (argparse::nargs_pattern::any)
      .remaining ();

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

  FailMode failMode = failModeFromOptions (options.get<std::string> ("-f"));

  std::string maxUlp = options.get<std::string> ("-m");

  gCheckExc = options.get<bool> ("-e");
  gCheckErrno = options.get<bool> ("-E");
  gComputeExc = gCheckExc || gCheckErrno;
  refimpls_compute_exc = gComputeExc ? 1 : 0;

  if (auto descFile = options.present ("-d"))
    handleDescription (*descFile, roundModes, failMode, maxUlp);
  else if (auto symbol = options.present ("-s"))
    {
      try
	{
	  handleList (*symbol,
		      options.get<std::vector<std::string> > ("values"),
		      roundModes, failMode, maxUlp);
	}
      catch (std::logic_error &e)
	{
	  error ("no values provided");
	}
    }
  else
    error ("no -d or -s provided");
}
