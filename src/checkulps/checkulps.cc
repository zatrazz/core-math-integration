//
// Copyright (c) Adhemerval Zanella. All rights reserved.
// Licensed under the MIT license. See LICENSE file in the project root for
// details.
//

#include <algorithm>
#include <iostream>
#include <limits>
#include <numbers>
#include <random>
#include <ranges>
#include <utility>

#include <argparse/argparse.hpp>

#include <cerrno>
#include <fenv.h>
#include <mpfr.h>
#include <omp.h>

#include "description.h"
#include "floatranges.h"
#include "floatsampler.h"
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

// Set by --summary: report only the maximum ULP found for each range/mode
// instead of the full histogram (or every checked value).
static bool gSummary = false;

// Print the ULP histogram body: every bin, or -- under --summary -- only the
// largest ULP distance seen (the last, highest key).
template <typename F>
static void
printUlpBins (const UlpAccumulator<F> &ulpacc, std::uint64_t ulptotal)
{
  auto bins = ulpacc.sorted ();
  if (gSummary)
    {
      if (bins.empty ())
	printlnTimestamp ("    max ulp: (no samples)");
      else
	{
	  const auto &last = *std::prev (bins.end ());
	  printlnTimestamp ("    max ulp {:g}  ({} of {}, {:6.2f}%)", last.first,
			    last.second, ulptotal,
			    ((double) last.second / (double) ulptotal) * 100.0);
	}
      return;
    }
  for (const auto &ulp : bins)
    printlnTimestamp ("    {:g}: {:16} {:6.2f}%", ulp.first, ulp.second,
		      ((double) ulp.second / (double) ulptotal) * 100.0);
}

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

  printUlpBins (ulpacc, ulptotal);
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

  printUlpBins (ulpacc, ulptotal);
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

  printUlpBins (ulpacc, ulptotal);
}

template <typename F>
static void
printAccumulator (const std::string_view &rndname,
		  const Description::FullRange &sample,
		  const UlpAccumulator<F> &ulpacc)
{
  const std::uint64_t ulptotal = ulpacc.total ();

  printlnTimestamp ("Checking rounding mode {:13}, {}", rndname, sample.name);

  printUlpBins (ulpacc, ulptotal);
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
// Distribution used to draw random inputs.  binade (uniform over representable
// floats) exercises every magnitude/binade, unlike real (uniform over the real
// interval), where a wide range's top binade dominates.
static floatsampler::Dist gDist = floatsampler::Dist::binade;
// True when either option needs the reference exception side-channel (errno
// expectations are derived from the expected exceptions).
static bool gComputeExc = false;
// Set by --special: force the special / corner-input check for every checked
// function, as if each description carried "special": true.
static bool gForceSpecial = false;
// Set by --sample: when >= 0, check only this (0-based) entry of the
// description's "samples" array instead of every sample.
static int gSampleIndex = -1;

// The canonical errno a correctly-rounded result with exceptions EXC and value
// EXPECTEDVALUE sets, following glibc: EDOM for a domain error (invalid) and
// ERANGE for a pole (divide-by-zero), for overflow when the result is infinite,
// and for underflow when the result flushes to zero.
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

// Whether GOT is an acceptable errno for a correctly-rounded result raising
// EXC with value EXPVAL.  EDOM (domain error) and ERANGE (pole) are required,
// but ISO C makes errno optional for a range error: on overflow or underflow a
// conforming libm may or may not set ERANGE (glibc's own choice varies with the
// function and whether the result reaches infinity / flushes to zero), so both
// 0 and ERANGE are accepted there to avoid spurious mismatches.
template <typename F>
static inline bool
errnoAcceptable (unsigned exc, F expVal, int got)
{
  if ((exc & (FE_OVERFLOW | FE_UNDERFLOW)) && !(exc & (FE_INVALID | FE_DIVBYZERO)))
    return got == 0 || got == ERANGE;
  return got == expectedErrno (exc, expVal);
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
  if (doErrno && !errnoAcceptable (expExc, expVal, gotErrno))
    reportErrnoMismatch (mk (), expectedErrno (expExc, expVal), gotErrno,
			 failmode);
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

  floatsampler::Sampler<FloatType> dist (gDist, sample.arg.start,
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

  floatsampler::Sampler<FloatType> dist (gDist, sample.arg.start,
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

  floatsampler::Sampler<FloatType> distX (gDist, sample.arg_x.start,
						  sample.arg_x.end);
  floatsampler::Sampler<FloatType> distY (gDist, sample.arg_y.start,
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

  floatsampler::Sampler<FloatType> distX (gDist, sample.arg_x.start,
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

      F maxUlp = -1;
      F maxInput = 0;
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
	  if (gSummary)
	    {
	      if (ret.ulp > maxUlp)
		{
		  maxUlp = ret.ulp;
		  maxInput = values[i];
		}
	    }
	  else if (!ret.checkFull ())
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
	  if (gCheckErrno
	      && !errnoAcceptable (expexc[i][idx], expected[i][idx], gotErrno))
	    reportErrnoMismatch (ret,
				 expectedErrno (expexc[i][idx], expected[i][idx]),
				 gotErrno, failmode);
	}
      if (gSummary && maxUlp >= 0)
	printlnTimestamp ("Checking rounding mode {:13}  max ulp {:1.0f}  "
			  "input={:#a}",
			  rnd.name, maxUlp, maxInput);
    }

  printlnTimestamp ("");
}

//
// Special / corner input checks.
//

template <typename F>
static std::vector<F>
specialValues ()
{
  return {
    F (0.0),
    -F (0.0),
    std::numeric_limits<F>::infinity (),
    -std::numeric_limits<F>::infinity (),
    std::numeric_limits<F>::quiet_NaN (),
    std::numeric_limits<F>::denorm_min (),
    -std::numeric_limits<F>::denorm_min (),
    std::numeric_limits<F>::min (),  // smallest positive normal
    -std::numeric_limits<F>::min (),
    std::numeric_limits<F>::max (),
    -std::numeric_limits<F>::max (),
    std::numbers::pi_v<F>, -std::numbers::pi_v<F>,
  };
}

static std::vector<long long int>
specialIntValues ()
{
  return { 0, 1, -1, 2, -2, 3, -3, 4, -4, 10, -10 };
}

// ---------------------------------------------------------------------------
// Worst cases for trigonometric range reduction.
//
// For a modulus C (pi/2, pi or pi/4) the hardest inputs for sin/cos/tan are
// the doubles x whose reduced argument x mod C is tiny: there the result is a
// small value obtained after massive cancellation, so a minute error in the
// reduction becomes a huge ULP error.  Such inputs form a measure-zero set and
// are never produced by uniform range sampling.  In the binade at exponent e a
// double is x = M * 2^(e-52) with M in [2^52, 2^53); the ones best
// approximating a multiple of C minimise ||M * g|| with g = frac(2^(e-52) / C),
// and those minimisers are exactly the convergents and semiconvergents of the
// continued fraction of g.  g is built with MPFR at high precision because the
// low bits of 1/C that survive the 2^(e-52) scaling drive the result (the
// Payne-Hanek phenomenon).
// ---------------------------------------------------------------------------

// For the candidate denominator M report the residual rho = |M*g - round(M*g)|
// (how close x = M*2^(e-52) is to a multiple of the modulus) and the parity of
// that nearest multiple.  Even parity means x is near k*C with k even, i.e. a
// zero of one function; odd parity is the shifted zero of the other -- both are
// worst cases, for different functions, so the caller keeps the worst of each.
// rho is returned as long double, whose 15-bit exponent preserves the ordering
// of residuals far below the double underflow threshold.
static void
reducedInfo (const mpfr_t g, unsigned long M, long double &rho, int &parity)
{
  mpfr_t t, r;
  mpfr_init2 (t, mpfr_get_prec (g));
  mpfr_init2 (r, mpfr_get_prec (g));
  mpfr_mul_ui (t, g, M, MPFR_RNDN);
  mpfr_round (r, t);
  parity = static_cast<int> (mpfr_get_uj (r, MPFR_RNDN) & 1u);
  mpfr_sub (t, t, r, MPFR_RNDN);
  mpfr_abs (t, t, MPFR_RNDN);
  rho = mpfr_get_ld (t, MPFR_RNDN);
  mpfr_clear (t);
  mpfr_clear (r);
}

static std::vector<double>
reductionWorstCases (const std::string &modulo, int exp_lo, int exp_hi,
		     uint64_t count)
{
  const uint64_t QMIN = 1ULL << 52;
  const uint64_t QMAX = (1ULL << 53) - 1;

  exp_lo = std::max (exp_lo, 1);
  exp_hi = std::min (exp_hi, 1023);
  if (count == 0)
    count = 1;

  // Precision must survive the 2^(e-52) scaling (which shifts up to ~971 bits
  // out) and still leave room for denominators up to 2^53 and the residual.
  const mpfr_prec_t prec
      = static_cast<mpfr_prec_t> (std::max (exp_hi, 64) + 256);

  // invC = 1 / C = k / pi, with k = 2 (pi/2), 1 (pi) or 4 (pi/4).
  unsigned long k = modulo == "pi/2" ? 2 : (modulo == "pi/4" ? 4 : 1);
  mpfr_t invC, pi, y, g;
  mpfr_init2 (invC, prec);
  mpfr_init2 (pi, prec);
  mpfr_init2 (y, prec);
  mpfr_init2 (g, prec);
  mpfr_const_pi (pi, MPFR_RNDN);
  mpfr_ui_div (invC, k, pi, MPFR_RNDN);

  std::vector<double> out;

  for (int e = exp_lo; e <= exp_hi; e++)
    {
      // g = frac(2^(e-52) * invC), the fractional part that governs ||M * g||.
      mpfr_mul_2si (y, invC, e - 52, MPFR_RNDN);
      mpfr_frac (g, y, MPFR_RNDN);
      if (mpfr_zero_p (g))
	continue;

      // The doubles closest to a multiple of C in the window are small integer
      // combinations of the continued-fraction convergent denominators q_k of
      // g.  In particular a q_k with tiny residual delta_k makes every small
      // multiple m*q_k have residual m*delta_k, still tiny -- so the worst
      // cases are not only the convergents/semiconvergents but the multiples of
      // convergents, offset by a few q_{k-1}.  Collect the convergent
      // denominators up to 2*QMAX, then enumerate those combinations that land
      // in the binade window.
      constexpr uint64_t kMultCap = 256; // largest multiplier per convergent
      const unsigned __int128 QLIM = (unsigned __int128) 2 * QMAX;

      std::vector<uint64_t> qs{ 1 }; // convergent denominators q_0, q_1, ...
      {
	uint64_t km2 = 0, km1 = 1; // q_{-1}, q_0
	mpfr_t t;
	mpfr_init2 (t, prec);
	mpfr_ui_div (t, 1, g, MPFR_RNDN); // t_1 = 1/g
	for (int step = 0; step < 4 * prec; step++)
	  {
	    if (mpfr_cmp_d (t, 9.0e18) >= 0) // a_k huge -> q_k overshoots window
	      break;
	    uint64_t ak = mpfr_get_uj (t, MPFR_RNDD);
	    unsigned __int128 qk = (unsigned __int128) ak * km1 + km2;
	    if (qk > QLIM)
	      break;
	    qs.push_back (static_cast<uint64_t> (qk));
	    km2 = km1;
	    km1 = static_cast<uint64_t> (qk);
	    mpfr_sub_ui (t, t, ak, MPFR_RNDN);
	    if (mpfr_zero_p (t))
	      break;
	    mpfr_ui_div (t, 1, t, MPFR_RNDN);
	  }
	mpfr_clear (t);
      }

      std::vector<uint64_t> cands{ QMIN, QMAX };
      static const int kOffsets[] = { 0, 1, -1, 2, -2 };
      for (std::size_t i = 1; i < qs.size (); i++)
	{
	  uint64_t qk = qs[i], qkm1 = qs[i - 1];
	  uint64_t mmax = QMAX / qk;
	  if (mmax > kMultCap)
	    continue; // tiny convergent: its multiples are covered by larger ones
	  for (uint64_t m = 1; m <= mmax; m++)
	    for (int o : kOffsets)
	      {
		__int128 d = static_cast<__int128> (m) * qk
			     + static_cast<__int128> (o) * qkm1;
		if (d >= static_cast<__int128> (QMIN)
		    && d <= static_cast<__int128> (QMAX))
		  cands.push_back (static_cast<uint64_t> (d));
	      }
	}

      std::sort (cands.begin (), cands.end ());
      cands.erase (std::unique (cands.begin (), cands.end ()), cands.end ());

      // Rank the window denominators by residual, then keep the smallest `count`
      // for each parity so that both the sin-type (even) and cos-type (odd)
      // worst cases of the binade are emitted.
      struct Cand
      {
	uint64_t M;
	long double rho;
	int parity;
      };
      std::vector<Cand> ranked;
      ranked.reserve (cands.size ());
      for (uint64_t M : cands)
	{
	  long double rho;
	  int parity;
	  reducedInfo (g, M, rho, parity);
	  ranked.push_back ({ M, rho, parity });
	}
      std::sort (ranked.begin (), ranked.end (),
		 [] (const Cand &a, const Cand &b) { return a.rho < b.rho; });

      uint64_t taken[2] = { 0, 0 };
      for (const Cand &c : ranked)
	{
	  if (taken[c.parity] >= count)
	    continue;
	  taken[c.parity]++;
	  double x = std::ldexp (static_cast<double> (c.M), e - 52);
	  // The function's worst input can be a near neighbour of the reduction
	  // worst case, so include a couple of ulps either side.
	  double xm2 = std::nextafter (std::nextafter (x, 0.0), 0.0);
	  double xm1 = std::nextafter (x, 0.0);
	  double xp1 = std::nextafter (x, 2.0 * x);
	  double xp2 = std::nextafter (xp1, 2.0 * x);
	  out.insert (out.end (), { xm2, xm1, x, xp1, xp2 });
	  if (taken[0] >= count && taken[1] >= count)
	    break;
	}
    }

  mpfr_clear (invC);
  mpfr_clear (pi);
  mpfr_clear (y);
  mpfr_clear (g);
  return out;
}

// Compare one already-evaluated result against the reference: report a value
// failure (honouring FailMode), then an exception and/or errno mismatch.
template <typename RET, typename F>
static void
reportListResult (const RET &ret, unsigned raised, unsigned expExc,
		  int gotErrno, F expVal, FailMode failmode)
{
  // Under --summary the caller reports the maximum ULP; suppress the per-value
  // lines here (exception/errno mismatches are still reported).
  if (gSummary)
    ; // no per-value output
  else if (!ret.checkFull ())
    {
      if (failmode == FailMode::FIRST || failmode == FailMode::ALL)
	{
	  printlnErrorTimestamp ("{}", ret);
	  if (failmode == FailMode::FIRST)
	    std::exit (EXIT_FAILURE);
	}
    }
  else
    printlnTimestamp ("{}", ret);

  if (gCheckExc && raised != expExc)
    reportExcMismatch (ret, expExc, raised, failmode);
  if (gCheckErrno && !errnoAcceptable (expExc, expVal, gotErrno))
    reportErrnoMismatch (ret, expectedErrno (expExc, expVal), gotErrno,
			 failmode);
}

// Two-argument explicit value-list check (used for the special cross product).
template <typename F>
static void
checkListFloatFloat (const std::vector<std::pair<F, F> > &values,
		     FuncFF<F> func, const FuncFFReference<F> &ref, F max_ulp,
		     const RoundSet &roundModes, FailMode failmode)
{
  const unsigned mask = maskFromRoundSet (roundModes);
  std::vector<std::array<F, REF_NRND> > expected (values.size ());
  std::vector<std::array<unsigned, REF_NRND> > expexc (values.size ());
  {
    RoundSetup<F> roundSetup (FE_TONEAREST);
    for (std::size_t i = 0; i < values.size (); i++)
      {
	ref (values[i].first, values[i].second, mask, expected[i].data ());
	if (gComputeExc)
	  captureExpExc (expexc[i].data ());
      }
  }

  for (const auto &rnd : roundModes)
    {
      int idx = refIndex (rnd.mode);
      RoundSetup<F> roundSetup (rnd.mode);
      F maxUlp = -1, maxX = 0, maxY = 0;
      for (std::size_t i = 0; i < values.size (); i++)
	{
	  if (gCheckExc)
	    feclearexcept (kDriverExcMask);
	  if (gCheckErrno)
	    errno = 0;
	  F computed = func (values[i].first, values[i].second);
	  unsigned raised
	      = gCheckExc ? (unsigned) fetestexcept (kDriverExcMask) : 0u;
	  int gotErrno = gCheckErrno ? errno : 0;
	  ResultFloatFloat<F> ret (rnd.mode, values[i].first, values[i].second,
				   computed, expected[i][idx], max_ulp);
	  if (gSummary && ret.ulp > maxUlp)
	    {
	      maxUlp = ret.ulp;
	      maxX = values[i].first;
	      maxY = values[i].second;
	    }
	  reportListResult (ret, raised, expexc[i][idx], gotErrno,
			    expected[i][idx], failmode);
	}
      if (gSummary && maxUlp >= 0)
	printlnTimestamp ("Checking rounding mode {:13}  max ulp {:1.0f}  "
			  "x={:#a} y={:#a}",
			  rnd.name, maxUlp, maxX, maxY);
    }
  printlnTimestamp ("");
}

// Float-and-integer explicit value-list check.
template <typename F>
static void
checkListFloatLLI (const std::vector<std::pair<F, long long int> > &values,
		   FuncFLLI<F> func, const FuncFLLIReference<F> &ref, F max_ulp,
		   const RoundSet &roundModes, FailMode failmode)
{
  const unsigned mask = maskFromRoundSet (roundModes);
  std::vector<std::array<F, REF_NRND> > expected (values.size ());
  std::vector<std::array<unsigned, REF_NRND> > expexc (values.size ());
  {
    RoundSetup<F> roundSetup (FE_TONEAREST);
    for (std::size_t i = 0; i < values.size (); i++)
      {
	ref (values[i].first, values[i].second, mask, expected[i].data ());
	if (gComputeExc)
	  captureExpExc (expexc[i].data ());
      }
  }

  for (const auto &rnd : roundModes)
    {
      int idx = refIndex (rnd.mode);
      RoundSetup<F> roundSetup (rnd.mode);
      F maxUlp = -1, maxX = 0;
      long long int maxN = 0;
      for (std::size_t i = 0; i < values.size (); i++)
	{
	  if (gCheckExc)
	    feclearexcept (kDriverExcMask);
	  if (gCheckErrno)
	    errno = 0;
	  F computed = func (values[i].first, values[i].second);
	  unsigned raised
	      = gCheckExc ? (unsigned) fetestexcept (kDriverExcMask) : 0u;
	  int gotErrno = gCheckErrno ? errno : 0;
	  ResultFloatLLI<F> ret (rnd.mode, values[i].first, values[i].second,
				 computed, expected[i][idx], max_ulp);
	  if (gSummary && ret.ulp > maxUlp)
	    {
	      maxUlp = ret.ulp;
	      maxX = values[i].first;
	      maxN = values[i].second;
	    }
	  reportListResult (ret, raised, expexc[i][idx], gotErrno,
			    expected[i][idx], failmode);
	}
      if (gSummary && maxUlp >= 0)
	printlnTimestamp ("Checking rounding mode {:13}  max ulp {:1.0f}  "
			  "x={:#a} n={}",
			  rnd.name, maxUlp, maxX, maxN);
    }
  printlnTimestamp ("");
}

// Single-input, two-output explicit value-list check (e.g. sincos).
template <typename F>
static void
checkListFloatpFloatp (const std::vector<F> &values, FuncFpFp<F> func,
		       const FuncFpFpReference<F> &ref, F max_ulp,
		       const RoundSet &roundModes, FailMode failmode)
{
  const unsigned mask = maskFromRoundSet (roundModes);
  std::vector<std::array<F, REF_NRND> > exp0 (values.size ()),
      exp1 (values.size ());
  std::vector<std::array<unsigned, REF_NRND> > expexc (values.size ());
  {
    RoundSetup<F> roundSetup (FE_TONEAREST);
    for (std::size_t i = 0; i < values.size (); i++)
      {
	ref (values[i], mask, exp0[i].data (), exp1[i].data ());
	if (gComputeExc)
	  captureExpExc (expexc[i].data ());
      }
  }

  for (const auto &rnd : roundModes)
    {
      int idx = refIndex (rnd.mode);
      RoundSetup<F> roundSetup (rnd.mode);
      F maxUlp = -1, maxInput = 0;
      for (std::size_t i = 0; i < values.size (); i++)
	{
	  if (gCheckExc)
	    feclearexcept (kDriverExcMask);
	  if (gCheckErrno)
	    errno = 0;
	  F computed0, computed1;
	  func (values[i], &computed0, &computed1);
	  unsigned raised
	      = gCheckExc ? (unsigned) fetestexcept (kDriverExcMask) : 0u;
	  int gotErrno = gCheckErrno ? errno : 0;
	  ResultFloatpFloatp<F> ret (rnd.mode, values[i], computed0, computed1,
				     exp0[i][idx], exp1[i][idx], max_ulp);
	  if (gSummary && ret.ulp > maxUlp)
	    {
	      maxUlp = ret.ulp;
	      maxInput = values[i];
	    }
	  reportListResult (ret, raised, expexc[i][idx], gotErrno,
			    exp0[i][idx], failmode);
	}
      if (gSummary && maxUlp >= 0)
	printlnTimestamp ("Checking rounding mode {:13}  max ulp {:1.0f}  "
			  "input={:#a}",
			  rnd.name, maxUlp, maxInput);
    }
  printlnTimestamp ("");
}

// Build the special cross product for a two-argument function.
template <typename F>
static std::vector<std::pair<F, F> >
specialPairs ()
{
  auto v = specialValues<F> ();
  std::vector<std::pair<F, F> > pairs;
  pairs.reserve (v.size () * v.size ());
  for (F x : v)
    for (F y : v)
      pairs.emplace_back (x, y);
  return pairs;
}

template <typename F>
static std::vector<std::pair<F, long long int> >
specialPairsLLI ()
{
  auto v = specialValues<F> ();
  auto n = specialIntValues ();
  std::vector<std::pair<F, long long int> > pairs;
  pairs.reserve (v.size () * n.size ());
  for (F x : v)
    for (long long int y : n)
      pairs.emplace_back (x, y);
  return pairs;
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
      else if (auto *psample
	       = std::get_if<Description::ReductionRange> (&sample))
	{
	  if constexpr (std::is_same_v<F, double>)
	    checkList (desc.FunctionName,
		       reductionWorstCases (psample->modulo, psample->exp_lo,
					    psample->exp_hi, psample->count),
		       func.first, func.second, max_ulp.value (), roundModes,
		       failmode);
	  else
	    error ("reduction sampling is only supported for double");
	}
      else
	error ("invalid sample type");
    }

  if (desc.CheckSpecial)
    checkList (desc.FunctionName, specialValues<F> (), func.first, func.second,
	       max_ulp.value (), roundModes, failmode);

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
      else if (auto *psample
	       = std::get_if<Description::ReductionRange> (&sample))
	{
	  if constexpr (std::is_same_v<F, double>)
	    checkListFloatpFloatp (
		reductionWorstCases (psample->modulo, psample->exp_lo,
				     psample->exp_hi, psample->count),
		func.first, func.second, max_ulp.value (), roundModes,
		failmode);
	  else
	    error ("reduction sampling is only supported for double");
	}
      else
	error ("invalid sample type");
    }

  if (desc.CheckSpecial)
    checkListFloatpFloatp (specialValues<F> (), func.first, func.second,
			   max_ulp.value (), roundModes, failmode);

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

  if (desc.CheckSpecial)
    checkListFloatFloat (specialPairs<F> (), func.first, func.second,
			 max_ulp.value (), roundModes, failmode);

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

  if (desc.CheckSpecial)
    checkListFloatLLI (specialPairsLLI<F> (), func.first, func.second,
		       max_ulp.value (), roundModes, failmode);

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

  if (gForceSpecial)
    desc.CheckSpecial = true;

  if (gSampleIndex >= 0)
    {
      if (static_cast<std::size_t> (gSampleIndex) >= desc.Samples.size ())
	error ("--sample {} out of range ({} has {} sample(s))", gSampleIndex,
	       desc.FunctionName, desc.Samples.size ());
      Description::SampleType chosen = desc.Samples[gSampleIndex];
      desc.Samples.assign (1, chosen);
    }

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

  options.add_argument ("--special", "-p")
      .help ("also check the special / corner inputs (signed zeros, "
	     "infinities, NaN, subnormal and normal extremes, domain edges, "
	     "and for multi-argument functions their cross product), as if the "
	     "description carried \"special\": true")
      .flag ();

  options.add_argument ("--distribution", "-D")
      .help ("random input distribution: 'binade' (uniform over representable "
	     "floats, exercises every magnitude) or 'real' (uniform over the "
	     "real interval)")
      .default_value (std::string ("binade"));

  options.add_argument ("--summary", "-S")
      .help ("report only the maximum ULP found for each range and rounding "
	     "mode, instead of the full histogram or every checked value")
      .flag ();

  options.add_argument ("--sample", "-i")
      .help ("check only this 0-based entry of the description's \"samples\" "
	     "array instead of every sample")
      .nargs (1);

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
  gForceSpecial = options.get<bool> ("-p");
  gSummary = options.get<bool> ("-S");

  if (auto s = options.present ("-i"))
    {
      try
	{
	  gSampleIndex = std::stoi (*s);
	}
      catch (const std::exception &)
	{
	  error ("invalid --sample index: {}", *s);
	}
      if (gSampleIndex < 0)
	error ("--sample index must be non-negative");
    }

  if (auto d = floatsampler::distFromString (options.get<std::string> ("-D")))
    gDist = *d;
  else
    error ("invalid distribution: {}", options.get<std::string> ("-D"));

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
