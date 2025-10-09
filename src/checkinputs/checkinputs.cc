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

template <typename F>
static void
check_f (const std::string &input, bool ignore_errors)
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

      if (auto n = floatrange::fromStr<F> (line); n.has_value ())
	(*it).numbers.push_back (n.value ());
      else if (ignore_errors)
	std::println ("line {} invalid number {}: {}", line_number, line,
		      n.error ());
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
      std::println ("{0:20}: min={1:a} ({1:g}) max={2:a} ({2:g}) count={3}",
		    w.name, *minmaxt.first, *minmaxt.second,
		    w.numbers.size ());
    }
}

template <typename F>
static void
check_f_f (const std::string &input, bool ignore_errors)
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

      if (auto n = floatrange::fromStr<F> (numbers[0]); n.has_value ())
	(*it).numbers_x.push_back (n.value ());
      else if (ignore_errors)
	std::println ("line {} invalid number {}: {}", line_number, numbers[0],
		      n.error ());
      else
	error ("line {} invalid number {}: {}", line_number, numbers[0],
	       n.error ());

      if (auto n = floatrange::fromStr<F> (numbers[1]); n.has_value ())
	(*it).numbers_y.push_back (n.value ());
      else if (ignore_errors)
	std::println ("line {} invalid number {}: {}", line_number, numbers[1],
		      n.error ());
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

  std::string type;
  options.add_argument ("--type", "-t")
      .help ("floating type to use")
      .default_value ("binary32")
      .store_into (type);

  bool bivariate;
  options.add_argument ("--bivariate", "-b")
      .help ("handle inputs as bivariate functions")
      .store_into (bivariate)
      .flag ();

  bool ignore_errors;
  options.add_argument ("--ignore_errors", "-i")
      .help ("Do not stop at first line parsing error")
      .store_into (ignore_errors)
      .flag ();

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
	check_f_f<float> (input, ignore_errors);
      else
	check_f<float> (input, ignore_errors);
    }
  else if (type == "binary64")
    {
      if (bivariate)
	check_f_f<double> (input, ignore_errors);
      else
	check_f<double> (input, ignore_errors);
    }
  else
    error ("invalid type {}", type);

  return 0;
}
