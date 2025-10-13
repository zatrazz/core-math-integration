//
// Copyright (c) Adhemerval Zanella. All rights reserved.
// Licensed under the MIT license. See LICENSE file in the project root for
// details.
//

#ifndef _DESCRIPTION_H
#define _DESCRIPTION_H

#include <cstdint>
#include <expected>
#include <string>
#include <vector>
#include <format>

#include "cxxcompat.h"

class Description
{
public:
  struct FullRange
  {
    std::string name;
    uint64_t start;
    uint64_t end;
  };

  template <typename F> struct ArgType
  {
    F start;
    F end;
  };

  template <typename F> struct Sample1Arg
  {
    ArgType<F> arg;
    uint64_t count;
  };

  template <typename F> struct Sample2Arg
  {
    ArgType<F> arg_x;
    ArgType<F> arg_y;
    uint64_t count;
  };

  template <typename F> struct Sample2ArgLli
  {
    ArgType<F> arg_x;
    ArgType<long long int> arg_y;
    uint64_t count;
  };

  std::expected<void, std::string> parse (const std::string &);

  // clang-format off
  typedef std::variant<Sample1Arg<float>,
		       Sample1Arg<double>,
		       Sample2Arg<float>,
                       Sample2Arg<double>,
		       Sample2ArgLli<float>,
                       Sample2ArgLli<double>,
		       FullRange>
      SampleType;
  // clang-format on

  std::string FunctionName;
  std::vector<SampleType> Samples;
};

template <>
struct std::formatter<Description::SampleType> : std::formatter<std::string>
{
  auto
  format (const Description::SampleType &s, std::format_context &ctx) const
  {
    return std::visit (
	[&] (auto &&arg) {
	  using T = std::decay_t<decltype (arg)>;
	  if constexpr (std::is_same_v<T, Description::Sample1Arg<float> >)
	    {
	      return std::format_to (ctx.out (), "Sample1Arg<float>: {}-{}",
				     arg.arg.start, arg.arg.end);
	    }
	  else if constexpr (std::is_same_v<T,
					    Description::Sample1Arg<double> >)
	    {
	      return std::format_to (ctx.out (), "Sample1Arg<double>: {}-{}",
				     arg.arg.start, arg.arg.end);
	    }
	  else if constexpr (std::is_same_v<T,
					    Description::Sample2Arg<float> >)
	    {
	      return std::format_to (ctx.out (),
				     "Sample2Arg<float>: {}-{} {}-{}",
				     arg.arg_x.start, arg.arg_x.end,
				     arg.arg_y.start, arg.arg_y.end);
	    }
	  else if constexpr (std::is_same_v<T,
					    Description::Sample2Arg<double> >)
	    {
	      return std::format_to (ctx.out (),
				     "Sample2Arg<double>: {}-{} {}-{}",
				     arg.arg_x.start, arg.arg_x.end,
				     arg.arg_y.start, arg.arg_y.end);
	    }
	  else if constexpr (std::is_same_v<
				 T, Description::Sample2ArgLli<float> >)
	    {
	      return std::format_to (ctx.out (),
				     "Sample2ArgLli<float>: {}-{} {}-{}",
				     arg.arg_x.start, arg.arg_x.end,
				     arg.arg_y.start, arg.arg_y.end);
	    }
	  else if constexpr (std::is_same_v<
				 T, Description::Sample2ArgLli<double> >)
	    {
	      return std::format_to (ctx.out (),
				     "Sample2ArgLli<double>: {}-{} {}-{}",
				     arg.arg_x.start, arg.arg_x.end,
				     arg.arg_y.start, arg.arg_y.end);
	    }
	  else if constexpr (std::is_same_v<T, Description::FullRange>)
	    {
	      return std::format_to (ctx.out (), "FullRange: {} {}-{}",
				     arg.name, arg.start, arg.end);
	    }
	},
	s);
  }
};

#endif
