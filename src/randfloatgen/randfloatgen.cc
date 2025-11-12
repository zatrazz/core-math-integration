//
// Copyright (c) Adhemerval Zanella. All rights reserved.
// Licensed under the MIT license. See LICENSE file in the project root for
// details.
//

#include <algorithm>
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
#include "wyhash64.h"

using namespace iohelper;
typedef wyhash64 rng_t;

static constexpr int kDefaultCount = 1000;

static rng_t
init_random_state (void)
{
  std::random_device rd{ "/dev/random" };
  return rng_t ((rng_t::state_type) rd () << 32 | rd ());
}

template <typename F>
static void
gen_f (const std::optional<std::string> &nameopt,
       const std::optional<std::string> &argsopt, F fstart, F fend, int count)
{
  std::string name;
  if (nameopt.has_value ())
    name = *nameopt;
  else
    name = "random";

  std::string args;
  if (argsopt.has_value ())
    std::println ("## args: {}", *argsopt);
  else
    {
      std::println ("## args: {}", floatrange::Limits<F>::name);
      std::println ("## ret: {}", floatrange::Limits<F>::name);
    }

  std::println ("## includes: math.h");
  std::println ("## name: workload-{}", name);
  std::println ("# Random inputs in [{:.2f},{:.2f}]", fstart, fend);

  rng_t rng = init_random_state ();
  std::uniform_real_distribution<F> d (fstart, fend);
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
	 const std::optional<std::string> &argsopt, F fstart0, F fend0,
	 F fstart1, F fend1, int count)
{
  std::string name;
  if (nameopt.has_value ())
    name = *nameopt;
  else
    name = "random";

  std::string args;
  if (argsopt.has_value ())
    args = *argsopt;
  else
    args = std::format ("## args {0}:{0}", floatrange::Limits<F>::name);

  std::println ("## ret: {}", floatrange::Limits<F>::name);
  std::println ("## includes: math.h");
  std::println ("## name: workload-{}", name);
  std::println (
      "# Random inputs with in in [{:.2f},{:.2f}] and y in [{:.2f},{:.2f}]",
      fstart0, fend0, fstart1, fend1);

  rng_t rng = init_random_state ();
  std::uniform_real_distribution<F> d0 (fstart0, fend0);
  std::uniform_real_distribution<F> d1 (fstart1, fend1);
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
	   const std::optional<std::string> &argsopt,
	   F fstart0, F fend0,
	   F fstart1, F fend1,
	   F fstart2, F fend2,
	   int count)
{
  std::string name;
  if (nameopt.has_value ())
    name = *nameopt;
  else
    name = "random";

  std::string args;
  if (argsopt.has_value ())
    std::println ("## args: {}", *argsopt);
  else
    std::println ("## args: {0}:{0}:{0}", floatrange::Limits<F>::name);

  std::println ("## ret: {}", floatrange::Limits<F>::name);
  std::println ("## includes: math.h");
  std::println ("## name: workload-{}", name);
  std::println (
      "# Random inputs with in in [{:.2f},{:.2f}], y in [{:.2f},{:.2f}], "
      "and z in [{:.2f},{:.2f}]",
      fstart0, fend0, fstart1, fend1, fstart2, fend2);

  rng_t rng = init_random_state ();
  std::uniform_real_distribution<F> d0 (fstart0, fend0);
  std::uniform_real_distribution<F> d1 (fstart1, fend1);
  std::uniform_real_distribution<F> d2 (fstart2, fend2);
  for (int i = 0; i < count; i++)
    {
      F f0 = d0 (rng);
      F f1 = d1 (rng);
      F f2 = d2 (rng);
      bool isnegative0 = f0 < 0.0;
      bool isnegative1 = f1 < 0.0;
      bool isnegative2 = f2 < 0.0;
      std::println ("{}0x{:a}, {}0x{:a}, {}0x{:a}",
		    isnegative0 ? "-" : "", std::fabs (f0),
		    isnegative1 ? "-" : "", std::fabs (f1),
		    isnegative2 ? "-" : "", std::fabs (f2));
    }
}

static void
handle_f (const std::optional<std::string> &nameopt,
	  const std::optional<std::string> &argsopt, const std::string &type,
	  int count, const std::vector<double> &ranges)
{
  if (type == "binary32")
    gen_f<float> (nameopt, argsopt, ranges[0], ranges[1], count);
  else if (type == "binary64")
    gen_f<double> (nameopt, argsopt, ranges[0], ranges[1], count);
  else
    error ("invalid type {}", type);
}

static void
handle_f_f (const std::optional<std::string> &nameopt,
	    const std::optional<std::string> &argsopt, const std::string &type,
	    int count, const std::vector<double> &ranges_x,
	    const std::vector<double> &ranges_y)
{
  if (type == "binary32")
    gen_f_f<float> (nameopt, argsopt, ranges_x[0], ranges_x[1], ranges_y[0],
		    ranges_y[1], count);
  else if (type == "binary64")
    gen_f_f<double> (nameopt, argsopt, ranges_x[0], ranges_x[1], ranges_y[0],
		     ranges_y[1], count);
  else
    error ("invalid type {}", type);
}

static void
handle_f_f_f (const std::optional<std::string> &nameopt,
	      const std::optional<std::string> &argsopt,
	      const std::string &type,
	      int count,
	      const std::vector<double> &ranges_x,
	      const std::vector<double> &ranges_y,
	      const std::vector<double> &ranges_z)
{
  if (type == "binary32")
    gen_f_f_f<float> (nameopt, argsopt, ranges_x[0], ranges_x[1], ranges_y[0],
		      ranges_y[1], ranges_z[0], ranges_z[1], count);
  else if (type == "binary64")
    gen_f_f_f<double> (nameopt, argsopt, ranges_x[0], ranges_x[1], ranges_y[0],
		     ranges_y[1], ranges_z[0], ranges_z[1],count);
  else
    error ("invalid type {}", type);
}

[[noreturn]] static inline void
error (const std::string &str)
{
  std::cerr << "error: " << str << '\n';
  std::exit (EXIT_FAILURE);
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

  // We always parse in binary64, even for binary32
  options.add_argument ("-x")
      .help ("range to use in the form '<start> <end>'")
      .required ()
      .nargs (2)
      .scan<'g', double> ();

  options.add_argument ("-y")
      .help ("range to use in the form '<start> <end>'")
      .nargs (2)
      .scan<'g', double> ();

  options.add_argument ("-z")
      .help ("range to use in the form '<start> <end>'")
      .nargs (2)
      .scan<'g', double> ();

  options.add_argument ("--count", "-c")
      .help (std::format ("numbers to generate (default {}", kDefaultCount))
      .default_value (kDefaultCount)
      .scan<'i', int>();

  options.add_argument ("--name", "-n");
  options.add_argument ("--args", "-a");

  try
    {
      options.parse_args (argc, argv);
    }
  catch (const std::runtime_error &err)
    {
      error ("{}", err.what ());
    }

  int count = options.get<int>("-c");

  std::optional<std::string> name;
  if (options.is_used ("--name"))
    name = options.get<std::string> ("--name");
  std::optional<std::string> args;
  if (options.is_used ("--args"))
    name = options.get<std::string> ("--args");

  auto range_x = options.get<std::vector<double> > ("-x");
  if (range_x[0] > range_x[1])
    error ("invalid range definitions [{},{}]", range_x[0], range_x[1]);

  if (!options.is_used ("-y"))
    handle_f (name, args, type, count, range_x);
  else
    {
      auto range_y = options.get<std::vector<double> > ("-y");
      if (range_y[0] > range_y[1])
	error ("invalid range definitions [{},{}]", range_y[0], range_y[1]);

      if (!options.is_used ("-z"))
	handle_f_f (name, args, type, count, range_x, range_y);
      else
	{
	  auto range_z = options.get<std::vector<double> >("-z");
	  if (range_z[0] > range_z[1])
	    error ("invalid range definitions [{},{}]", range_z[0], range_z[1]);

	  handle_f_f_f (name, args, type, count, range_x, range_y, range_z);
	}
    }
}
