//
// Copyright (c) Adhemerval Zanella. All rights reserved.
// Licensed under the MIT license. See LICENSE file in the project root for
// details.
//

#ifndef UINT128_H
#define UINT128_H

#include "config.h"
#include <cstdint>

#if HAVE_UINT128_T
typedef __uint128_t uint128_t;
#else
// Limited 128-bit integer support when __uint128_t is not nativelly
// supported (used on wyhash64).

class uint128_t
{
  std::uint64_t low{ 0 };
  std::uint64_t high{ 0 };

  constexpr
  uint128_t (std::uint64_t h, std::uint64_t l)
  {
    low = l;
    high = h;
  }

  static constexpr std::uint64_t TOPBIT = UINT64_C (1) << 63;

public:
  constexpr
  uint128_t ()
  {
  }

  constexpr
  uint128_t (int n)
  {
    low = static_cast<std::uint64_t> (n);
    high = -static_cast<std::uint64_t> (n < 0);
  }

  constexpr
  uint128_t (unsigned n)
      : low{ n }
  {
  }
  constexpr
  uint128_t (unsigned long n)
      : low{ n }
  {
  }
  constexpr
  uint128_t (unsigned long long n)
      : low{ n }
  {
  }

  constexpr uint128_t (const uint128_t &) = default;
  constexpr uint128_t (uint128_t &&) = default;
  constexpr uint128_t &operator= (const uint128_t &) = default;
  constexpr uint128_t &operator= (uint128_t &&) = default;

  constexpr explicit
  operator std::uint64_t () const
  {
    return low;
  }

  constexpr bool
  operator> (uint128_t that) const
  {
    return that < *this;
  }

  constexpr bool
  operator<= (uint128_t that) const
  {
    return !(*this > that);
  }

  constexpr bool
  operator>= (uint128_t that) const
  {
    return that <= *this;
  }

  constexpr bool
  operator< (uint128_t that) const
  {
    return high < that.high || (high == that.high && low < that.low);
  }

  constexpr bool
  operator== (uint128_t that) const
  {
    return low == that.low && high == that.high;
  }

  constexpr uint128_t
  operator^ (uint128_t that) const
  {
    return { high ^ that.high, low ^ that.low };
  }

  constexpr uint128_t
  operator<< (uint128_t that) const
  {
    if (that >= 128)
      return {};
    else if (that == 0)
      return *this;
    else
      {
	std::uint64_t n{ that.low };
	if (n >= 64)
	  return { low << (n - 64), 0 };
	else
	  return { (high << n) | (low >> (64 - n)), low << n };
      }
  }

  constexpr uint128_t
  operator>> (uint128_t that) const
  {
    if (that >= 128)
      return {};
    else if (that == 0)
      return *this;
    else
      {
	std::uint64_t n{ that.low };
	if (n >= 64)
	  return { 0, high >> (n - 64) };
	else
	  return { high >> n, (high << (64 - n)) | (low >> n) };
      }
  }

  constexpr uint128_t
  operator+ (uint128_t that) const
  {
    std::uint64_t lower{ (low & ~TOPBIT) + (that.low & ~TOPBIT) };
    bool carry{ ((lower >> 63) + (low >> 63) + (that.low >> 63)) > 1 };
    return { high + that.high + carry, low + that.low };
  }

  constexpr uint128_t
  operator* (uint128_t that) const
  {
    std::uint64_t mask32{ 0xffffffff };
    if (high == 0 && that.high == 0)
      {
	std::uint64_t x0{ low & mask32 }, x1{ low >> 32 };
	std::uint64_t y0{ that.low & mask32 }, y1{ that.low >> 32 };
	uint128_t x0y0{ x0 * y0 }, x0y1{ x0 * y1 };
	uint128_t x1y0{ x1 * y0 }, x1y1{ x1 * y1 };
	return x0y0 + ((x0y1 + x1y0) << 32) + (x1y1 << 64);
      }
    else
      {
	std::uint64_t x0{ low & mask32 }, x1{ low >> 32 }, x2{ high & mask32 },
	    x3{ high >> 32 };
	std::uint64_t y0{ that.low & mask32 }, y1{ that.low >> 32 },
	    y2{ that.high & mask32 }, y3{ that.high >> 32 };
	uint128_t x0y0{ x0 * y0 }, x0y1{ x0 * y1 }, x0y2{ x0 * y2 },
	    x0y3{ x0 * y3 };
	uint128_t x1y0{ x1 * y0 }, x1y1{ x1 * y1 }, x1y2{ x1 * y2 };
	uint128_t x2y0{ x2 * y0 }, x2y1{ x2 * y1 };
	uint128_t x3y0{ x3 * y0 };
	return x0y0 + ((x0y1 + x1y0) << 32) + ((x0y2 + x1y1 + x2y0) << 64)
	       + ((x0y3 + x1y2 + x2y1 + x3y0) << 96);
      }
  }
};
#endif

#endif
