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
      sample_type_t;
  // clang-format on

  std::string FunctionName;
  std::vector<sample_type_t> Samples;
};

#endif
