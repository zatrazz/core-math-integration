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

class description_t
{
public:
  struct full_t
  {
    std::string name;
    uint64_t start;
    uint64_t end;
  };

  template <typename F> struct arg_t
  {
    F start;
    F end;
  };

  template <typename F> struct sample_f
  {
    arg_t<F> arg;
    uint64_t count;
  };

  template <typename F> struct sample_f_f
  {
    arg_t<F> arg_x;
    arg_t<F> arg_y;
    uint64_t count;
  };

  template <typename F> struct sample_f_lli
  {
    arg_t<F> arg_x;
    arg_t<long long int> arg_y;
    uint64_t count;
  };

  std::expected<void, std::string> parse (const std::string &);

  // clang-format off
  typedef std::variant<sample_f<float>,
		       sample_f<double>,
		       sample_f_f<float>,
                       sample_f_f<double>,
		       sample_f_lli<float>,
                       sample_f_lli<double>,
		       full_t>
      sample_type_t;
  // clang-format on

  std::string function;
  std::vector<sample_type_t> samples;
};

#endif
