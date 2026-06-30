//
// Copyright (c) Adhemerval Zanella. All rights reserved.
// Licensed under the MIT license. See LICENSE file in the project root for
// details.
//

#include <algorithm>
#include <bit>
#include <climits>
#include <format>
#include <iostream>
#include <random>
#include <ranges>
#include <string_view>
#include <vector>

#include <argparse/argparse.hpp>

#include "floatranges.h"
#include "iohelper.h"
#include "strhelper.h"
#include "wyhash64.h"

using namespace iohelper;
typedef wyhash64 rng_t;

static constexpr int kDefaultCount = 1000;

// Sampling distribution for the generated values.
//
//  - real: uniform over the real interval [a,b] (std::uniform_real_-
//	        distribution).  Samples are spread by real-number measure, so the
//	        top binade dominates and small magnitudes/subnormals almost never
//	        appear.  Good for "typical input" benchmark workloads.
//
//  - binade: uniform over the representable floats in [a,b].  Each ULP is
//	          equally likely, so every binade (exponent) gets the same number of
//	          samples, naturally exercising small values and subnormals.  Good
//	          for accuracy/ULP coverage.  Also fully reproducible across
//	          toolchains, unlike uniform_real_distribution.
enum class Dist
{
  real,
  binade,
};

template <typename F> struct FloatBits;
template <> struct FloatBits<float>
{
  using type = uint32_t;
};
template <> struct FloatBits<double>
{
  using type = uint64_t;
};

// Map an IEEE float to a monotonically increasing unsigned key, giving a total
// order across the sign bit (most negative -> 0, ... -0, +0 ..., +max ->
// all-ones).  Sampling a uniform integer between two keys therefore samples
// uniformly over the representable values of the range.
template <typename F>
static typename FloatBits<F>::type
to_ordered (F f)
{
  using U = typename FloatBits<F>::type;
  constexpr U sign = U (1) << (sizeof (U) * CHAR_BIT - 1);
  U b = std::bit_cast<U> (f);
  return (b & sign) ? ~b : (b | sign);
}

template <typename F>
static F
from_ordered (typename FloatBits<F>::type o)
{
  using U = typename FloatBits<F>::type;
  constexpr U sign = U (1) << (sizeof (U) * CHAR_BIT - 1);
  U b = (o & sign) ? (o & ~sign) : ~o;
  return std::bit_cast<F> (b);
}

template <typename F> class Sampler
{
  Dist dist;
  std::uniform_real_distribution<F> real;
  std::uniform_int_distribution<typename FloatBits<F>::type> bits;

public:
  Sampler (Dist d, F a, F b)
      : dist (d), real (a, b), bits (to_ordered<F> (a), to_ordered<F> (b))
  {
  }

  template <typename RNG>
  F
  operator() (RNG &rng)
  {
    return dist == Dist::real ? real (rng) : from_ordered<F> (bits (rng));
  }
};

static rng_t
init_random_state (void)
{
  std::random_device rd{ "/dev/random" };
  return rng_t ((rng_t::state_type) rd () << 32 | rd ());
}

template <typename F>
static void
gen_f (const std::optional<std::string> &nameopt,
       const std::optional<std::string> &argsopt,
       F fstart, F fend, int count,
       bool append, Dist dist)
{
  std::string name;
  if (nameopt.has_value ())
    name = *nameopt;
  else
    name = "random";

  if (!append)
    {
      std::string args;
      if (argsopt.has_value ())
	std::println ("## args: {}", *argsopt);
      else
	{
	  std::println ("## args: {}", floatrange::Limits<F>::name);
	  std::println ("## ret: {}", floatrange::Limits<F>::name);
	}
      std::println ("## includes: math.h");
    }

  std::println ("## name: workload-{}", name);
  std::println ("# Random inputs in [{:.2f},{:.2f}]", fstart, fend);

  rng_t rng = init_random_state ();
  Sampler<F> d (dist, fstart, fend);
  for (int i = 0; i < count; i++)
    {
      F f = d (rng);
      bool isnegative = f < 0.0;
      std::println ("{}0x{:a}", isnegative ? "-" : "", std::fabs (f));
    }
}

template <typename F>
static void
gen_f_f (const std::optional<std::string> &nameopt,
	 const std::optional<std::string> &argsopt,
	 F fstart0, F fend0, F fstart1, F fend1, int count,
	 bool append, Dist dist)
{
  std::string name;
  if (nameopt.has_value ())
    name = *nameopt;
  else
    name = "random";

  if (!append)
    {
      if (argsopt.has_value ())
	std::println ("## args: {}", *argsopt);
      else
	std::println ("## args: {0}:{0}", floatrange::Limits<F>::name);
      std::println ("## ret: {}", floatrange::Limits<F>::name);
      std::println ("## includes: math.h");
    }
  std::println ("## name: workload-{}", name);
  std::println (
      "# Random inputs with x in [{:.2f},{:.2f}] and y in [{:.2f},{:.2f}]",
      fstart0, fend0, fstart1, fend1);

  rng_t rng = init_random_state ();
  Sampler<F> d0 (dist, fstart0, fend0);
  Sampler<F> d1 (dist, fstart1, fend1);
  for (int i = 0; i < count; i++)
    {
      F f0 = d0 (rng);
      F f1 = d1 (rng);
      bool isnegative0 = f0 < 0.0;
      bool isnegative1 = f1 < 0.0;
      std::println ("{}0x{:a}, {}0x{:a}", isnegative0 ? "-" : "",
		    std::fabs (f0), isnegative1 ? "-" : "", std::fabs (f1));
    }
}

template <typename F>
static void
gen_f_f_f (const std::optional<std::string> &nameopt,
	   const std::optional<std::string> &argsopt, F fstart0, F fend0,
	   F fstart1, F fend1, F fstart2, F fend2, int count,
	   bool append, Dist dist)
{
  std::string name;
  if (nameopt.has_value ())
    name = *nameopt;
  else
    name = "random";

  if (!append)
    {
      std::string args;
      if (argsopt.has_value ())
	std::println ("## args: {}", *argsopt);
      else
	std::println ("## args: {0}:{0}:{0}", floatrange::Limits<F>::name);

      std::println ("## ret: {}", floatrange::Limits<F>::name);
      std::println ("## includes: math.h");
    }
  std::println ("## name: workload-{}", name);
  std::println (
      "# Random inputs with in in [{:.2f},{:.2f}], y in [{:.2f},{:.2f}], "
      "and z in [{:.2f},{:.2f}]",
      fstart0, fend0, fstart1, fend1, fstart2, fend2);

  rng_t rng = init_random_state ();
  Sampler<F> d0 (dist, fstart0, fend0);
  Sampler<F> d1 (dist, fstart1, fend1);
  Sampler<F> d2 (dist, fstart2, fend2);
  for (int i = 0; i < count; i++)
    {
      F f0 = d0 (rng);
      F f1 = d1 (rng);
      F f2 = d2 (rng);
      bool isnegative0 = f0 < 0.0;
      bool isnegative1 = f1 < 0.0;
      bool isnegative2 = f2 < 0.0;
      std::println ("{}0x{:a}, {}0x{:a}, {}0x{:a}", isnegative0 ? "-" : "",
		    std::fabs (f0), isnegative1 ? "-" : "", std::fabs (f1),
		    isnegative2 ? "-" : "", std::fabs (f2));
    }
}

[[noreturn]] static inline void
error (const std::string &str)
{
  std::cerr << "error: " << str << '\n';
  std::exit (EXIT_FAILURE);
}

static inline std::string
adjustSignal (const std::string &s)
{
  return s.starts_with ("\\-") ? "-" + s.substr (2) : s;
}

template <typename F>
static const std::vector<F>
rangeStrToFloat (const std::vector<std::string> &values)
{
  auto floatview = values | std::views::transform ([] (const std::string &s) {
		     auto trimmed = adjustSignal (std::string(strhelper::trim (s)));
		     if (trimmed == "-pi")
		       return -std::numbers::pi_v<F>;
		     else if (trimmed == "pi")
		       return std::numbers::pi_v<F>;
		     else if (trimmed == "2pi")
		       return F(2.0) * std::numbers::pi_v<F>;
		     else if (trimmed == "m2pi")
		       return -F(2.0) * std::numbers::pi_v<F>;
		     else if (trimmed == "min")
		       return std::numeric_limits<F>::min ();
		     else if (trimmed == "mmin")
		       return -std::numeric_limits<F>::min ();
		     else if (trimmed == "max")
		       return std::numeric_limits<F>::max ();
		     else if (trimmed == "mmax")
		       return -std::numeric_limits<F>::max ();
		     if (auto n = floatrange::fromStr<F> (s); n.has_value ())
		       return n.value ();
		     else
		       error ("invalid number: {}", s);
		   });

  return std::vector<F> (floatview.begin (), floatview.end ());
}

template <typename F>
static void
handleType (const argparse::ArgumentParser &options,
	    const std::optional<std::string> &name,
	    const std::optional<std::string> &args,
	    int count,
	    bool append, Dist dist)
{
  auto range_x
      = rangeStrToFloat<F> (options.get<std::vector<std::string> > ("-x"));
  if (range_x[0] > range_x[1])
    error ("invalid range definitions [{},{}]", range_x[0], range_x[1]);

  if (!options.is_used ("-y"))
    gen_f<F> (name, args, range_x[0], range_x[1], count, append, dist);
  else
    {
      auto range_y
	  = rangeStrToFloat<F> (options.get<std::vector<std::string> > ("-y"));
      if (range_y[0] > range_y[1])
	error ("invalid range definitions [{},{}]", range_y[0], range_y[1]);

      if (!options.is_used ("-z"))
	gen_f_f<F> (name, args, range_x[0], range_x[1], range_y[0], range_y[1],
		    count, append, dist);
      else
	{
	  auto range_z = rangeStrToFloat<F> (
	      options.get<std::vector<std::string> > ("-z"));
	  if (range_z[0] > range_z[1])
	    error ("invalid range definitions [{},{}]", range_z[0],
		   range_z[1]);

	  gen_f_f_f<F> (name, args, range_x[0], range_x[1], range_y[0],
			range_y[1], range_z[0], range_z[1], count, append,
			dist);
	}
    }
}

int
main (int argc, char *argv[])
{
  argparse::ArgumentParser options ("randfloatgen");

  std::string type;
  options.add_argument ("--type", "-t")
      .help (std::format ("floating type to use"))
      .default_value ("binary32")
      .store_into (type);

  options.add_argument ("-x")
      .help ("range to use in the form '<start> <end>'")
      .required ()
      .nargs (2);

  options.add_argument ("-y")
      .help ("range to use in the form '<start> <end>'")
      .nargs (2);

  options.add_argument ("-z")
      .help ("range to use in the form '<start> <end>'")
      .nargs (2);

  options.add_argument ("--count", "-c")
      .help (std::format ("numbers to generate (default {}", kDefaultCount))
      .default_value (kDefaultCount)
      .scan<'i', int> ();

  options.add_argument ("--dist", "-d")
      .help ("sampling distribution: 'real' (uniform over the real interval) "
	     "or 'binade' (uniform over representable floats)")
      .default_value (std::string ("real"));

  options.add_argument ("--name", "-n");
  options.add_argument ("--args", "-a");

  options.add_argument ("--append")
     .default_value (false)
     .implicit_value (true);

  try
    {
      options.parse_args (argc, argv);
    }
  catch (const std::runtime_error &err)
    {
      error ("{}", err.what ());
    }

  int count = options.get<int> ("-c");

  std::optional<std::string> name;
  if (options.is_used ("--name"))
    name = options.get<std::string> ("--name");
  std::optional<std::string> args;
  if (options.is_used ("--args"))
    args = options.get<std::string> ("--args");

  bool append = options.get<bool>("--append");

  Dist dist;
  std::string distStr = options.get<std::string> ("--dist");
  if (distStr == "real")
    dist = Dist::real;
  else if (distStr == "binade")
    dist = Dist::binade;
  else
    error ("invalid distribution: {}", distStr);

  if (type == "binary32")
    handleType<float> (options, name, args, count, append, dist);
  else if (type == "binary64")
    handleType<double> (options, name, args, count, append, dist);
  else
    error ("invalid type{}", type);
}
