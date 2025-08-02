#include <fstream>
#include <numbers>
#include <ranges>

#include <nlohmann/json.hpp>

#include "cxxcompat.h"
#include "description.h"
#include "floatranges.h"
#include "refimpls.h"
#include "strhelper.h"

static const std::array full_ranges_options = { "normal", "subnormal" };

static std::vector<std::string>
split_with_ranges (const std::string_view &s, std::string_view delimiter)
{
  std::vector<std::string> tokens;
  for (const auto &subrange : s | std::views::split (delimiter))
    tokens.push_back (std::string (subrange.begin (), subrange.end ()));
  return tokens;
}

#define TRY(expr)                                                             \
  ({                                                                          \
    auto __result = (expr);                                                   \
    if (!__result.has_value ())                                               \
      return std::unexpected (__result.error ());                             \
    __result.value ();                                                        \
  })

static std::expected<description_t::full_t, std::string>
handle_full (refimpls::func_type_t functype, const std::string &name)
{
  if (name == "normal")
    {
      switch (functype)
        {
        case refimpls::func_type_t::f32_f:
          return description_t::full_t{
            "positive normal", float_ranges_t::limits<float>::plus_normal_min,
            float_ranges_t::limits<float>::plus_normal_max
          };
        case refimpls::func_type_t::f64_f:
          return description_t::full_t{
            "negative normal", float_ranges_t::limits<double>::plus_normal_min,
            float_ranges_t::limits<double>::plus_normal_max
          };
        default:
          std::unreachable ();
        }
    }
  else if (name == "subnormal")
    {
      switch (functype)
        {
        case refimpls::func_type_t::f32_f:
          return description_t::full_t{
            "positive subnormal",
            float_ranges_t::limits<float>::plus_subnormal_min,
            float_ranges_t::limits<float>::plus_subnormal_max
          };
        case refimpls::func_type_t::f64_f:
          return description_t::full_t{
            "negative subnormal",
            float_ranges_t::limits<double>::plus_subnormal_min,
            float_ranges_t::limits<double>::plus_subnormal_max
          };
        default:
          std::unreachable ();
        }
    }
  return std::unexpected (std::format ("invalid full range: {}", name));
}

template <typename F>
inline std::expected<F, std::string>
parse_range (const std::string &str)
{
  std::string_view trimmed = strhelper::trim (str);
  if (trimmed == "-pi")
    return -std::numbers::pi_v<F>;
  else if (trimmed == "pi")
    return std::numbers::pi_v<F>;
  else if (trimmed == "2pi")
    return 2.0 * std::numbers::pi_v<F>;
  else if (trimmed == "min")
    return std::numeric_limits<F>::min ();
  else if (trimmed == "-min")
    return -std::numeric_limits<F>::min ();
  else if (trimmed == "max")
    return std::numeric_limits<F>::max ();
  else if (trimmed == "-max")
    return -std::numeric_limits<F>::max ();
  return float_ranges_t::from_str<F> (std::string (trimmed));
}

template <>
inline std::expected<long long int, std::string>
parse_range (const std::string &str)
{
  long long int r;
  auto result = std::from_chars (str.data (), str.data () + str.size (), r);
  if (result.ec == std::errc{})
    return r;
  else if (result.ec == std::errc::invalid_argument)
    return std::unexpected (std::format ("invalid number: {}", str));
  else if (result.ec == std::errc::result_out_of_range)
    return std::unexpected (
        std::format ("number out of range for long long: {}", str));
  std::unreachable ();
}

static std::expected<description_t::sample_type_t, std::string>
handle_1_arg (refimpls::func_type_t functype, const std::string &start,
              const std::string &end, uint64_t count)
{
  switch (functype)
    {
    case refimpls::func_type_t::f32_f:
      {

        return description_t::sample_f<float>{
          TRY (parse_range<float> (start)), TRY (parse_range<float> (end)),
          count
        };
      }
    case refimpls::func_type_t::f64_f:
      return description_t::sample_f<double>{
        TRY (parse_range<double> (start)), TRY (parse_range<double> (end)),
        count
      };
    default:
      std::unreachable ();
    }
}

static std::expected<description_t::sample_type_t, std::string>
handle_2_arg (refimpls::func_type_t functype, const std::string &start_x,
              const std::string &end_x, const std::string &start_y,
              const std::string &end_y, uint64_t count)
{
  switch (functype)
    {
    case refimpls::func_type_t::f32_f_f:
      return description_t::sample_f_f<float>{
        TRY (parse_range<float> (start_x)), TRY (parse_range<float> (end_x)),
        TRY (parse_range<float> (start_y)), TRY (parse_range<float> (end_y)),
        count
      };
    case refimpls::func_type_t::f64_f_f:
      return description_t::sample_f_f<double>{
        TRY (parse_range<double> (start_x)), TRY (parse_range<double> (end_x)),
        TRY (parse_range<double> (start_y)), TRY (parse_range<double> (end_y)),
        count
      };
    case refimpls::func_type_t::f32_f_lli:
      return description_t::sample_f_lli<float>{
        TRY (parse_range<float> (start_x)), TRY (parse_range<float> (end_x)),
        TRY (parse_range<long long int> (start_y)),
        TRY (parse_range<long long int> (end_y)), count
      };
    case refimpls::func_type_t::f64_f_lli:
      return description_t::sample_f_lli<double>{
        TRY (parse_range<double> (start_x)), TRY (parse_range<double> (end_x)),
        TRY (parse_range<long long int> (start_y)),
        TRY (parse_range<long long int> (end_y)), count
      };
    default:
      std::unreachable ();
    }
}

template<class... Ts>
struct overloaded : Ts... { using Ts::operator()...; };


std::expected<void, std::string>
description_t::parse (const std::string &fname)
{
  std::ifstream file (fname);

  nlohmann::json data;
  try
   {
     data = nlohmann::json::parse (file);
   }
  catch (nlohmann::json::parse_error& ex)
   {
     return std::unexpected (std::format ("error parsing file {} at position {}",
                             fname, ex.byte));
   }

  this->function = data["function"];
  auto functype = refimpls::get_func_type (function);
  if (!functype)
    return std::unexpected (std::format ("invalid function: {}", function));

  if (data.contains ("full"))
    {
      auto fulldata
          = split_with_ranges (data["full"].get<std::string_view> (), ",");
      for (const auto &f : fulldata)
        {
          auto full = handle_full (functype.value (), f);
          if (!full)
            return std::unexpected (full.error ());
          this->samples.push_back (full.value ());
        }
    }
  else if (data.contains ("samples"))
    {
      for (const auto &r : data["samples"])
        {
          if (r.contains ("x") && r.contains ("y"))
            {
              if (r["x"].size () != 2 || r["y"].size () != 2)
                return std::unexpected (
                    std::format ("invalid sample size: ({}, {}) (expected 2)",
                                 r["x"].size (), r["y"].size ()));

              auto sample = TRY (handle_2_arg (
                  functype.value (), r["x"][0].template get<std::string> (),
                  r["x"][1].template get<std::string> (),
                  r["y"][0].template get<std::string> (),
                  r["y"][1].template get<std::string> (),
                  r["count"].get<uint64_t> ()));
              this->samples.push_back (sample);
            }
          else if (r.contains ("x"))
            {
              if (r["x"].size () != 2)
                return std::unexpected (std::format (
                    "invalid sample size: {} (expected 2)", r["x"].size ()));

              auto sample = TRY (handle_1_arg (
                  functype.value (), r["x"][0].template get<std::string> (),
                  r["x"][1].template get<std::string> (),
                  r["count"].get<uint64_t> ()));

              this->samples.push_back (sample);
            }
          else
            return std::unexpected (
                std::format ("invalid sample definition {}", r.dump ()));
        }
    }
  else
    return std::unexpected ("no samples found");

  return {};
}
