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
	  auto fields = strhelper::splitWithRanges (line.substr (2), ":");
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
	  auto fields = strhelper::splitWithRanges (line.substr (2), ":");
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

      auto numbers = strhelper::splitWithRanges (line, ",");
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
      std::println (
	  "x = [{0:a}, {1:a}] ({0:g}, {1:g})",
	  *minmaxt_x.first, *minmaxt_x.second);

      auto minmaxt_y
	  = std::minmax_element (w.numbers_y.begin (), w.numbers_y.end ());
      std::println (
	  "y = [{0:a}, {1:a}] ({0:g}, {1:g})",
	  *minmaxt_y.first, *minmaxt_y.second);
    }
}

template <typename F>
static void
check_f_f_f (const std::string &input, bool ignore_errors)
{
  std::ifstream file (input);
  if (!file.is_open ())
    error ("opening file {}", input);

  struct workload_t
  {
    std::string name;
    std::vector<F> numbers_x;
    std::vector<F> numbers_y;
    std::vector<F> numbers_z;
  };

  std::vector<workload_t> workloads = { { "default" } };
  auto it = workloads.begin ();

  std::string line;
  for (int line_number = 1; std::getline (file, line); line_number++)
    {
      if (line.starts_with ("##"))
	{
	  auto fields = strhelper::splitWithRanges (line.substr (2), ":");
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

      auto numbers = strhelper::splitWithRanges (line, ",");
      if (numbers.size () != 3)
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

      if (auto n = floatrange::fromStr<F> (numbers[2]); n.has_value ())
	(*it).numbers_z.push_back (n.value ());
      else if (ignore_errors)
	std::println ("line {} invalid number {}: {}", line_number, numbers[2],
		      n.error ());
      else
	error ("line {} invalid number {}: {}", line_number, numbers[2],
	       n.error ());
    }
  if (workloads.empty ())
    return;

  for (const auto &w : workloads)
    {
      if (w.numbers_x.size () == 0 || w.numbers_y.size () == 0
	  || w.numbers_z.size( ) == 0)
	continue;

      auto minmaxt_x
	  = std::minmax_element (w.numbers_x.begin (), w.numbers_x.end ());
      std::println (
	  "x = [{0:a}, {1:a}] ({0:g}, {1:g})",
	  *minmaxt_x.first, *minmaxt_x.second);

      auto minmaxt_y
	  = std::minmax_element (w.numbers_y.begin (), w.numbers_y.end ());
      std::println (
	  "y = [{0:a}, {1:a}] ({0:g}, {1:g})",
	  *minmaxt_y.first, *minmaxt_y.second);

      auto minmaxt_z
	  = std::minmax_element (w.numbers_z.begin (), w.numbers_z.end ());
      std::println (
	  "z = [{0:a}, {1:a}] ({0:g}, {1:g})",
	  *minmaxt_z.first, *minmaxt_z.second);
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

  options.add_argument ("--args", "-n")
      .help ("number of arguments")
      .default_value (1)
      .scan<'i', int>()
      .required();

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

  int nargs = options.get<int>("-n");
  if (nargs < 0 || nargs > 3)
    error ("invalid number of arguments ({})", nargs);

  if (type == "binary32")
    {
      switch (nargs)
	{
	case 1: check_f<float> (input, ignore_errors); break;
	case 2: check_f_f<float> (input, ignore_errors); break;
	case 3: check_f_f_f<float> (input, ignore_errors); break;
	}
    }
  else if (type == "binary64")
    {
      switch (nargs)
	{
	case 1: check_f<double> (input, ignore_errors); break;
	case 2: check_f_f<double> (input, ignore_errors); break;
	case 3: check_f_f_f<double> (input, ignore_errors); break;
	}
    }
  else if (type == "binary96")
    {
      switch (nargs)
	{
	case 1: check_f<long double> (input, ignore_errors); break;
	case 2: check_f_f<long double> (input, ignore_errors); break;
	case 3: check_f_f_f<long double> (input, ignore_errors); break;
	}
    }
  else
    error ("invalid type {}", type);

  return 0;
}
