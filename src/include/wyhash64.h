//
// Copyright (c) Adhemerval Zanella. All rights reserved.
// Licensed under the MIT license. See LICENSE file in the project root for
// details.
//

#ifndef WYHASH64_H
#define WYHASH64_H

#include <cstdint>
#include <limits>
#include "uint128.h"

// From https://github.com/lemire/testingRNG/blob/master/source/wyhash.h with
// a C++ wrapper to be used with random distribution classes.

class wyhash64
{
public:
  typedef std::uint64_t result_type;
  typedef std::uint64_t state_type;

  static constexpr size_t state_size = sizeof (std::uint64_t);
  static constexpr size_t default_seed = 0;

  explicit wyhash64 () { seed (default_seed); }

  explicit wyhash64 (result_type __sd) { seed (__sd); }

  void
  seed (result_type __sd = default_seed)
  {
    state = __sd;
  }

  static constexpr result_type
  min ()
  {
    return std::numeric_limits<result_type>::min ();
  }

  static constexpr result_type
  max ()
  {
    return std::numeric_limits<result_type>::max ();
  }

  result_type
  operator() ()
  {
    state += 0x60bee2bee120fc15ull;
    uint128_t tmp = (uint128_t)state * 0xa3b195354a39b70dull;
    std::uint64_t m1 = std::uint64_t ((tmp >> 64) ^ tmp);
    tmp = (uint128_t)m1 * 0x1b03738712fad5c9ull;
    std::uint64_t m2 = std::uint64_t ((tmp >> 64) ^ tmp);
    return m2;
  }

private:
  std::uint64_t state;
};

#endif
