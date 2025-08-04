#include <algorithm>
#include <climits>
#include <format>
#include <iostream>
#include <random>
#include <ranges>
#include <string_view>
#include <vector>

#include <argparse/argparse.hpp>

#include "cxxcompat.h"
#include "wyhash64.h"

typedef wyhash64 rng_t;

template <typename T> struct ftypeinfo
{
  static std::string name ();
  static T fromstring (const std::string &);
};

template <> struct ftypeinfo<float>
{
  using type = float;
  static std::string
  name ()
  {
    return "float";
  }
  static float
  fromstring (const std::string &s)
  {
    return std::stof (s);
  }
};

template <> struct ftypeinfo<double>
{
  using type = double;
  static std::string
  name ()
  {
    return "double";
  }
  static double
  fromstring (const std::string &s)
  {
    return std::stod (s);
  }
};

static rng_t
init_random_state (void)
{
  std::random_device rd{ "/dev/random" };
  return rng_t ((rng_t::state_type)rd () << 32 | rd ());
}

template <typename ftypeinfo>
static void
gen_range_binary (const std::optional<std::string> &nameopt,
                  const std::optional<std::string> &argsopt,
                  const std::string &start, const std::string &end, int count)
{
  using ftype = typename ftypeinfo::type;

  ftype fstart = ftypeinfo::fromstring (start);
  ftype fend = ftypeinfo::fromstring (end);

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
      std::println ("{0}:{0}", ftypeinfo::name ());
      std::println ("## ret: {}", ftypeinfo::name ());
    }

  std::println ("## includes: math.h");
  std::println ("## name: workload-{}", name);
  std::println ("# Random inputs in [{:.2f},{:.2f}]", fstart, fend);

  rng_t rng = init_random_state ();
  std::uniform_real_distribution<ftype> d (fstart, fend);
  for (int i = 0; i < count; i++)
    {
      ftype f = d (rng);
      bool isnegative = f < 0.0;
      std::println ("{}0x{:a}", isnegative ? "-" : "", std::fabs (f));
    }
}

template <typename ftypeinfo>
static void
gen_range_binary_bivariate (const std::optional<std::string> &nameopt,
                            const std::optional<std::string> &argsopt,
                            const std::string &start0, const std::string &end0,
                            const std::string &start1, const std::string &end1,
                            int count)
{
  using ftype = typename ftypeinfo::type;

  ftype fstart0 = ftypeinfo::fromstring (start0);
  ftype fend0 = ftypeinfo::fromstring (end0);
  ftype fstart1 = ftypeinfo::fromstring (start1);
  ftype fend1 = ftypeinfo::fromstring (end1);

  std::string name;
  if (nameopt.has_value ())
    name = *nameopt;
  else
    name = "random";

  std::string args;
  if (argsopt.has_value ())
    args = *argsopt;
  else
    args = std::format ("{0}:{0}", ftypeinfo::name ());

  std::println ("## args: {0}:{0}", ftypeinfo::name ());
  std::println ("## ret: {}", ftypeinfo::name ());
  std::println ("## includes: math.h");
  std::println ("## name: workload-{}", name);
  std::println (
      "# Random inputs with in in [{:.2f},{:.2f}] and y in [{:.2f},{:.2f}]",
      fstart0, fend0, fstart1, fend1);

  rng_t rng = init_random_state ();
  std::uniform_real_distribution<ftype> d0 (fstart0, fend0);
  std::uniform_real_distribution<ftype> d1 (fstart1, fend1);
  for (int i = 0; i < count; i++)
    {
      ftype f0 = d0 (rng);
      ftype f1 = d1 (rng);
      bool isnegative0 = f0 < 0.0;
      bool isnegative1 = f1 < 0.0;
      std::println ("{}0x{:a}, {}0x{:a}", isnegative0 ? "-" : "",
                    std::fabs (f0), isnegative1 ? "-" : "", std::fabs (f1));
    }
}

static std::vector<std::string>
splitWithRanges (const std::string &s, std::string_view delimiter)
{
  std::vector<std::string> tokens;
  for (const auto &subrange : s | std::views::split (delimiter))
    tokens.push_back (std::string (subrange.begin (), subrange.end ()));
  return tokens;
}

static int
handle_univariate (const std::optional<std::string> &nameopt,
                   const std::optional<std::string> &argsopt,
                   const std::string &type, int count, std::string &range)
{
  std::vector<std::string> rangeFields = splitWithRanges (range, ":");
  if (rangeFields.size () != 2)
    {
      std::cout << "error: invalid range (" << range << ")\n";
      return 1;
    }

  if (type == "binary32")
    gen_range_binary<ftypeinfo<float> > (nameopt, argsopt, rangeFields[0],
                                         rangeFields[1], count);
  else if (type == "binary64")
    gen_range_binary<ftypeinfo<double> > (nameopt, argsopt, rangeFields[0],
                                          rangeFields[1], count);

  return 0;
}

static int
handle_bivariate (const std::optional<std::string> &nameopt,
                  const std::optional<std::string> &argsopt,
                  const std::string &type, int count, std::string &range)
{
  std::vector<std::string> ranges = splitWithRanges (range, " ");
  if (ranges.size () != 2)
    {
      std::cout << "error: invalid range (" << range << ")\n";
      return 1;
    }

  std::vector<std::string> range0 = splitWithRanges (ranges[0], ":");
  std::vector<std::string> range1 = splitWithRanges (ranges[1], ":");
  if (range0.size () != 2 || range1.size () != 2)
    {
      std::cout << "error: invalid range (" << range << ")\n";
      return 1;
    }

  if (type == "binary32")
    gen_range_binary_bivariate<ftypeinfo<float> > (
        nameopt, argsopt, range0[0], range0[1], range1[0], range1[1], count);
  else if (type == "binary64")
    gen_range_binary_bivariate<ftypeinfo<double> > (
        nameopt, argsopt, range0[0], range0[1], range1[0], range1[1], count);

  return 0;
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

  std::string type = "binary32";
  options.add_argument ("--type", "-t")
      .help ("floating type to use")
      .required ()
      .store_into (type);

  std::string range = "0.0:10.0";
  options.add_argument ("--range", "-r")
      .help ("range to use in the form '<start>:<end>'")
      .required ()
      .store_into (range);

  int count = 1000;
  options.add_argument ("--count", "-c")
      .help ("numbers to generate")
      .required ()
      .store_into (count);

  bool bivariate = false;
  options.add_argument ("--bivariate", "-b")
      .help ("handle inputs as bivariate functions")
      .store_into (bivariate);

  options.add_argument ("--name", "-n");
  options.add_argument ("--args", "-a");

  try
    {
      options.parse_args (argc, argv);
    }
  catch (const std::runtime_error &err)
    {
      error (std::string (err.what ()));
    }

  std::optional<std::string> name;
  if (options.is_used ("--name"))
    name = options.get<std::string> ("--name");
  std::optional<std::string> args;
  if (options.is_used ("--args"))
    name = options.get<std::string> ("--args");

  if (!bivariate)
    return handle_univariate (name, args, type, count, range);
  else
    return handle_bivariate (name, args, type, count, range);
}
