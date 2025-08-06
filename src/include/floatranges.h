//
// Copyright (c) Adhemerval Zanella. All rights reserved.
// Licensed under the MIT license. See LICENSE file in the project root for
// details.
//

#ifndef _FLOATRANGES_H
#define _FLOATRANGES_H

#include <cstdint>
#include <charconv>
#include <expected>
#include <format>
#include <string>

#include "cxxcompat.h"

namespace float_ranges_t
{

template <typename F, F (*conv) (const std::string &, std::size_t *)>
inline std::expected<F, std::string>
__from_str (const std::string &sv)
{
  try
    {
      size_t pos;
      F r = conv (sv, &pos);
      if (pos == sv.length ())
        return r;

      return std::unexpected (
          std::format ("invalid float conversion: {}", sv));
    }
  catch (const std::invalid_argument &e)
    {
      return std::unexpected (
          std::format ("invalid float conversion: {}", e.what ()));
    }
  catch (const std::out_of_range &e)
    {
      return std::unexpected (
          std::format ("number out of range: {}", e.what ()));
    }
}

template <typename T>
inline std::expected<T, std::string> from_str (const std::string &sv);

template <>
inline std::expected<float, std::string>
from_str (const std::string &sv)
{
  return __from_str<float, std::stof> (sv);
}

template <>
inline std::expected<double, std::string>
from_str (const std::string &sv)
{
  return __from_str<double, std::stod> (sv);
}

// Information class used to generate full ranges, mainly for testing
// all binary32 normal and subnormal numbers.
template <typename T> struct limits
{
};

template <> struct limits<float>
{
  static constexpr uint64_t plus_normal_min = 0x00800000;
  static constexpr uint64_t plus_normal_max = 0x7F7FFFFF;
  static constexpr uint64_t plus_subnormal_min = 0x00000001;
  static constexpr uint64_t plus_subnormal_max = 0x007FFFFF;
  static constexpr uint64_t neg_normal_min = 0x80800000;
  static constexpr uint64_t neg_normal_max = 0xFF7FFFFF;
  static constexpr uint64_t neg_subnormal_min = 0x80000001;
  static constexpr uint64_t neg_subnormal_max = 0x807FFFFF;

  static constexpr float
  from_integer (uint32_t u)
  {
    union
    {
      uint32_t u;
      float f;
    } r = { .u = u };
    return r.f;
  }
};

template <> struct limits<double>
{
  static constexpr uint64_t plus_normal_min = UINT64_C (0x0008000000000000);
  static constexpr uint64_t plus_normal_max = UINT64_C (0x7FEFFFFFFFFFFFFF);
  static constexpr uint64_t plus_subnormal_min = UINT64_C (0x0000000000000001);
  static constexpr uint64_t plus_subnormal_max = UINT64_C (0x000FFFFFFFFFFFFF);
  static constexpr uint64_t neg_normal_min = UINT64_C (0x8008000000000000);
  static constexpr uint64_t neg_normal_max = UINT64_C (0xFFEFFFFFFFFFFFFF);
  static constexpr uint64_t neg_subnormal_min = UINT64_C (0x8000000000000001);
  static constexpr uint64_t neg_subnormal_max = UINT64_C (0x800FFFFFFFFFFFFF);

  static constexpr double
  from_integer (uint64_t u)
  {
    union
    {
      uint64_t u;
      double f;
    } r = { .u = u };
    return r.f;
  }
};

} // namespace float_ranges_t

#endif
