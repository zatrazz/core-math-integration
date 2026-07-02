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

namespace floatrange
{

template <typename F, F (*conv) (const std::string &, std::size_t *)>
inline std::expected<F, std::string>
__fromStr (const std::string &sv)
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
inline std::expected<T, std::string> fromStr (const std::string &sv);

static inline std::string
removeTrailing (const std::string &s, char trailing)
{
  if (!s.empty ()
      && (s.back () == trailing || s.back () == std::toupper (trailing)))
    return s.substr (0, s.size () - 1);
  return s;
}

constexpr char ascii_tolower(char c) noexcept {
  return (c >= 'A' && c <= 'Z') ? c - 'A' + 'a' : c;
}

// Check if S is an IEEE infinity in the form [+]inf (case insensitive).
static inline bool
isInfinityToken (std::string_view v)
{
 if (!v.empty() && (v.front() == '+' || v.front() == '-'))
   v.remove_prefix(1);

  auto ieq = [](std::string_view a, std::string_view b) noexcept {
    // b assumed lowercase ASCII
    return std::ranges::equal(a, b, [](char x, char y) {
      return ascii_tolower(x) == y;
    });
  };
  return ieq (v, "inf") || ieq (v, "infinity");
}

template <>
inline std::expected<float, std::string>
fromStr (const std::string &sv)
{
  std::string rt = sv;
  if (rt.starts_with ("\\-"))
    rt = "-" + rt.substr (2);
  if (!isInfinityToken (rt))
    rt = removeTrailing (rt, 'f');
  return __fromStr<float, std::stof> (rt);
}

template <>
inline std::expected<double, std::string>
fromStr (const std::string &sv)
{
  if (sv.starts_with ("\\-"))
    return __fromStr<double, std::stod> ("-" + sv.substr (2));
  return __fromStr<double, std::stod> (sv);
}

template <>
inline std::expected<long double, std::string>
fromStr (const std::string &sv)
{
  std::string rt = sv;
  if (rt.starts_with ("\\-"))
    rt = "-" + rt.substr (2);
  if (!isInfinityToken (rt))
    rt = removeTrailing (rt, 'l');
  return __fromStr<long double, std::stold> (rt);
}

// Information class used to generate full ranges, mainly for testing
// all binary32 normal and subnormal numbers.
template <typename T> struct Limits
{
};

template <> struct Limits<float>
{
  static constexpr uint64_t PlusNormalMin = 0x00800000;
  static constexpr uint64_t PlusNormalMax = 0x7F7FFFFF;
  static constexpr uint64_t PlusSubnormalMin = 0x00000001;
  static constexpr uint64_t PlusSubnormalMax = 0x007FFFFF;
  static constexpr uint64_t NegNormalMin = 0x80800000;
  static constexpr uint64_t NegNormalMax = 0xFF7FFFFF;
  static constexpr uint64_t NegSubnormalMin = 0x80000001;
  static constexpr uint64_t NegSubnormalMax = 0x807FFFFF;

  static constexpr float
  from (uint32_t u)
  {
    union
    {
      uint32_t u;
      float f;
    } r = { .u = u };
    return r.f;
  }

  static constexpr const char name[] = "float";
};

template <> struct Limits<double>
{
  static constexpr uint64_t PlusNormalMin = UINT64_C (0x0010000000000000);
  static constexpr uint64_t PlusNormalMax = UINT64_C (0x7FEFFFFFFFFFFFFF);
  static constexpr uint64_t PlusSubnormalMin = UINT64_C (0x0000000000000001);
  static constexpr uint64_t PlusSubnormalMax = UINT64_C (0x000FFFFFFFFFFFFF);
  static constexpr uint64_t NegNormalMin = UINT64_C (0x8010000000000000);
  static constexpr uint64_t NegNormalMax = UINT64_C (0xFFEFFFFFFFFFFFFF);
  static constexpr uint64_t NegSubnormalMin = UINT64_C (0x8000000000000001);
  static constexpr uint64_t NegSubnormalMax = UINT64_C (0x800FFFFFFFFFFFFF);

  static constexpr double
  from (uint64_t u)
  {
    union
    {
      uint64_t u;
      double f;
    } r = { .u = u };
    return r.f;
  }

  static constexpr const char name[] = "double";
};

} // namespace float_ranges_t

#endif
