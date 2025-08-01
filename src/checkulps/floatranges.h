#ifndef _FLOATRANGES_H
#define _FLOATRANGES_H

#include <charconv>
#include <expected>
#include <string>

namespace float_ranges_t
{

template <typename T>
static std::expected<T, std::string>
from_str (const std::string_view &sv)
{
  T ret;
  auto [ptr, ec] = std::from_chars (sv.data (), sv.data () + sv.size (), ret);
  if (ec == std::errc{})
    return ret;
  return std::unexpected (std::format ("invalid float conversion: {}", sv));
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
