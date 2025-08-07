//
// Copyright (c) Adhemerval Zanella. All rights reserved.
// Licensed under the MIT license. See LICENSE file in the project root for
// details.
//

#include <algorithm>
#include <fstream>
#include <limits>
#include <string>

#include "floatranges.h"
#include "iohelper.h"
#include "strhelper.h"

#include <argparse/argparse.hpp>

using namespace iohelper;

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

template <> struct ftypeinfo<long double>
{
  using type = long double;
  static std::string
  name ()
  {
    return "long double";
  }
  static long double
  fromstring (const std::string &s)
  {
    return std::stold (s);
  }
};

template <typename F>
static void
check_f (const std::string &input)
{
  std::ifstream file (input);
  if (!file.is_open ())
    error ("opening file {}", input);

  struct workload_t
  {
    std::string name;
    std::vector<F> numbers;
  };

  std::vector<workload_t> workloads = { { "default" } };
  auto it = workloads.begin ();

  std::string line;
  for (int line_number = 1; std::getline (file, line); line_number++)
    {
      if (line.starts_with ("##"))
        {
          auto fields = strhelper::split_with_ranges (line.substr (2), ":");
          if (fields.size () != 2)
            error ("line {} invalid directive: {}", line_number, line);

          if (strhelper::trim (fields[0]).starts_with ("name"))
            {
              workloads.push_back (workload_t{ strhelper::trim (fields[1]) });
              it = workloads.end () - 1;
            }

          continue;
        }

      // Skip blank lines and comments.
      strhelper::trim (line);
      if (line.empty () || line.starts_with ("#"))
        continue;

      if (auto n = float_ranges_t::from_str<F> (line); n.has_value ())
        (*it).numbers.push_back (n.value ());
      else
        error ("line {} invalid number {}: {}", line_number, line, n.error ());
    }

  if (workloads.empty ())
    return;

  for (const auto &w : workloads)
    {
      if (w.numbers.size () == 0)
        continue;
      auto minmaxt
          = std::minmax_element (w.numbers.begin (), w.numbers.end ());
      std::println ("{0:20}: min={1:a} ({1:g})  max={2:a} ({2:g})", w.name,
                    *minmaxt.first, *minmaxt.second);
    }
}

template <typename F>
static void
check_f_f (const std::string &input)
{
  std::ifstream file (input);
  if (!file.is_open ())
    error ("opening file {}", input);

  struct workload_t
  {
    std::string name;
    std::vector<F> numbers_x;
    std::vector<F> numbers_y;
  };

  std::vector<workload_t> workloads = { { "default" } };
  auto it = workloads.begin ();

  std::string line;
  for (int line_number = 1; std::getline (file, line); line_number++)
    {
      if (line.starts_with ("##"))
        {
          auto fields = strhelper::split_with_ranges (line.substr (2), ":");
          if (fields.size () < 2)
            error ("line {} invalid directive: {}", line_number, line);

          if (strhelper::trim (fields[0]).starts_with ("name"))
            {
              workloads.push_back (workload_t{ strhelper::trim (fields[1]) });
              it = workloads.end () - 1;
            }

          continue;
        }

      // Skip blank lines and comments.
      strhelper::trim (line);
      if (line.empty () || line.starts_with ("#"))
        continue;

      auto numbers = strhelper::split_with_ranges (line, ",");
      if (numbers.size () != 2)
        error ("line {} not enough numbers: {}", line_number, line);

      if (auto n = float_ranges_t::from_str<F> (numbers[0]); n.has_value ())
        (*it).numbers_x.push_back (n.value ());
      else
        error ("line {} invalid number {}: {}", line_number, numbers[0],
               n.error ());

      if (auto n = float_ranges_t::from_str<F> (numbers[1]); n.has_value ())
        (*it).numbers_y.push_back (n.value ());
      else
        error ("line {} invalid number {}: {}", line_number, numbers[1],
               n.error ());
    }
  if (workloads.empty ())
    return;

  for (const auto &w : workloads)
    {
      if (w.numbers_x.size () == 0 || w.numbers_y.size () == 0)
        continue;
      auto minmaxt_x
          = std::minmax_element (w.numbers_x.begin (), w.numbers_x.end ());
      auto minmaxt_y
          = std::minmax_element (w.numbers_y.begin (), w.numbers_y.end ());

      std::println (
          "min={0:a}, {1:a} ({0:g}, {1:g})  max={2:a}, {3:a} ({2:g}, {3:g})",
          *minmaxt_x.first, *minmaxt_x.second, *minmaxt_y.first,
          *minmaxt_y.second);
    }
}

int
main (int argc, char *argv[])
{
  argparse::ArgumentParser options ("checkinputs");

  std::string type = "binary32";
  options.add_argument ("--type", "-t")
      .help ("floating type to use")
      .store_into (type);

  bool bivariate = false;
  options.add_argument ("--bivariate", "-b")
      .help ("handle inputs as bivariate functions")
      .store_into (bivariate);

  std::string input;
  options.add_argument ("input")
      .help ("glibc benchtest input file to parse")
      .nargs (1)
      .store_into (input)
      .required ();

  try
    {
      options.parse_args (argc, argv);
    }
  catch (const std::runtime_error &err)
    {
      error (std::string (err.what ()));
    }

  if (type == "binary32")
    {
      if (bivariate)
        check_f_f<float> (input);
      else
        check_f<float> (input);
    }
  else if (type == "binary64")
    {
      if (bivariate)
        check_f_f<double> (input);
      else
        check_f<double> (input);
    }

  return 0;
}
