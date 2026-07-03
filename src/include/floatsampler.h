//
// Copyright (c) Adhemerval Zanella. All rights reserved.
// Licensed under the MIT license. See LICENSE file in the project root for
// details.
//

#ifndef FLOATSAMPLER_H
#define FLOATSAMPLER_H

#include <algorithm>
#include <bit>
#include <climits>
#include <cstdint>
#include <optional>
#include <random>
#include <string_view>

namespace floatsampler
{

// Sampling distribution for random floating-point inputs.
//
//  - real: uniform over the real interval [a,b] (std::uniform_real_-
//	        distribution).  Samples are spread by real-number measure, so the
//	        top binade dominates and small magnitudes/subnormals almost never
//	        appear.  Good for "typical input" benchmark workloads.
//
//  - binade: uniform over the representable floats in [a,b].  Each ULP is
//	          equally likely, so every binade (exponent) gets the same number of
//	          samples, naturally exercising small values and subnormals.  Good
//	          for accuracy/ULP coverage.  Also fully reproducible across
//	          toolchains, unlike uniform_real_distribution.
enum class Dist
{
  real,
  binade,
};

template <typename F> struct FloatBits;
template <> struct FloatBits<float>
{
  using type = std::uint32_t;
};
template <> struct FloatBits<double>
{
  using type = std::uint64_t;
};

// Map an IEEE float to a monotonically increasing unsigned key, giving a total
// order across the sign bit (most negative -> 0, ... -0, +0 ..., +max ->
// all-ones).  Sampling a uniform integer between two keys therefore samples
// uniformly over the representable values of the range.
template <typename F>
typename FloatBits<F>::type
to_ordered (F f)
{
  using U = typename FloatBits<F>::type;
  constexpr U sign = U (1) << (sizeof (U) * CHAR_BIT - 1);
  U b = std::bit_cast<U> (f);
  return (b & sign) ? ~b : (b | sign);
}

template <typename F>
F
from_ordered (typename FloatBits<F>::type o)
{
  using U = typename FloatBits<F>::type;
  constexpr U sign = U (1) << (sizeof (U) * CHAR_BIT - 1);
  U b = (o & sign) ? (o & ~sign) : ~o;
  return std::bit_cast<F> (b);
}

template <typename F> class Sampler
{
  Dist dist;
  std::uniform_real_distribution<F> real;
  std::uniform_int_distribution<typename FloatBits<F>::type> bits;

public:
  Sampler (Dist d, F a, F b)
      : dist (d), real (a, b),
	bits (std::min (to_ordered<F> (a), to_ordered<F> (b)),
	      std::max (to_ordered<F> (a), to_ordered<F> (b)))
  {
  }

  template <typename RNG>
  F
  operator() (RNG &rng)
  {
    return dist == Dist::real ? real (rng) : from_ordered<F> (bits (rng));
  }
};

// Parse a distribution name ("real" or "binade"); std::nullopt if invalid.
inline std::optional<Dist>
distFromString (std::string_view s)
{
  if (s == "real")
    return Dist::real;
  if (s == "binade")
    return Dist::binade;
  return std::nullopt;
}

} // namespace floatsampler

#endif
