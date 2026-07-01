//
// Copyright (c) Adhemerval Zanella. All rights reserved.
// Licensed under the MIT license. See LICENSE file in the project root for
// details.
//

#include <algorithm>
#include <fstream>
#include <limits>
#include <string>
#include <vector>

#include "floatranges.h"
#include "iohelper.h"
#include "strhelper.h"

#include <argparse/argparse.hpp>

using namespace iohelper;

// Parse a glibc benchtest input file and report, for each workload, the
// min/max range of every argument column.
//
// The number of columns is taken from the command line (nargs) but is
// overridden by a "## args: T:U:..." directive when present, so a file can
// be inspected without knowing its arity up front.
template <typename F>
static void
check (const std::string &input, int nargs, bool ignore_errors)
{
  std::ifstream file (input);
  if (!file.is_open ())
    error ("opening file {}", input);

  struct workload_t
  {
    std::string name;
    std::vector<std::vector<F>> columns;
  };

  std::vector<workload_t> workloads;
  int ncols = nargs;

  // Return the workload the current data lines belong to, creating an
  // unnamed "default" one on demand for files without a "## name:"
  // directive.
  auto current = [&] () -> workload_t & {
    if (workloads.empty ())
      workloads.push_back ({ "default", std::vector<std::vector<F>> (ncols) });
    return workloads.back ();
  };

  std::string line;
  for (int line_number = 1; std::getline (file, line); line_number++)
    {
      if (line.starts_with ("##"))
	{
	  auto fields = strhelper::splitWithRanges (line.substr (2), ":");
	  if (fields.size () < 2)
	    error ("line {} invalid directive: {}", line_number, line);

	  auto key = strhelper::trim (fields[0]);
	  if (key.starts_with ("name"))
	    workloads.push_back ({ std::string (strhelper::trim (fields[1])),
				   std::vector<std::vector<F>> (ncols) });
	  else if (key.starts_with ("args"))
	    // "## args: T:U:..." lists one type per argument.
	    ncols = fields.size () - 1;

	  continue;
	}

      // Skip blank lines and comments.
      strhelper::trim (line);
      if (line.empty () || line.starts_with ("#"))
	continue;

      auto numbers = strhelper::splitWithRanges (line, ",");
      if (static_cast<int> (numbers.size ()) != ncols)
	{
	  if (ignore_errors)
	    {
	      std::println ("line {} expected {} number(s), found {}: {}",
			    line_number, ncols, numbers.size (), line);
	      continue;
	    }
	  error ("line {} expected {} number(s), found {}: {}", line_number,
		 ncols, numbers.size (), line);
	}

      workload_t &w = current ();
      if (static_cast<int> (w.columns.size ()) < ncols)
	w.columns.resize (ncols);

      for (int i = 0; i < ncols; i++)
	{
	  auto n = floatrange::fromStr<F> (
	      std::string (strhelper::trim (numbers[i])));
	  if (n.has_value ())
	    w.columns[i].push_back (n.value ());
	  else if (ignore_errors)
	    std::println ("line {} invalid number {}: {}", line_number,
			  numbers[i], n.error ());
	  else
	    error ("line {} invalid number {}: {}", line_number, numbers[i],
		   n.error ());
	}
    }

  static constexpr const char *labels[] = { "x", "y", "z" };
  for (const auto &w : workloads)
    {
      if (w.columns.empty () || w.columns[0].empty ())
	continue;

      std::println ("{:20}: count={}", w.name, w.columns[0].size ());
      for (size_t i = 0; i < w.columns.size (); i++)
	{
	  if (w.columns[i].empty ())
	    continue;
	  auto mm
	      = std::minmax_element (w.columns[i].begin (), w.columns[i].end ());
	  const char *label = i < std::size (labels) ? labels[i] : "arg";
	  std::println ("  {0} = [{1:a}, {2:a}] ({1:g}, {2:g})", label,
			*mm.first, *mm.second);
	}
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
      .help ("number of arguments (overridden by a '## args:' directive)")
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
  if (nargs < 1 || nargs > 3)
    error ("invalid number of arguments ({})", nargs);

  if (type == "binary32")
    check<float> (input, nargs, ignore_errors);
  else if (type == "binary64")
    check<double> (input, nargs, ignore_errors);
  else if (type == "binary96")
    check<long double> (input, nargs, ignore_errors);
  else
    error ("invalid type {}", type);

  return 0;
}
