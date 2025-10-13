//
// Copyright (c) Adhemerval Zanella. All rights reserved.
// Licensed under the MIT license. See LICENSE file in the project root for
// details.
//

#include <fstream>
#include <numbers>
#include <ranges>

#include <nlohmann/json.hpp>

#include "cxxcompat.h"
#include "description.h"
#include "floatranges.h"
#include "refimpls.h"
#include "strhelper.h"

static std::vector<std::string>
SplitWithRanges (const std::string_view &s, std::string_view delimiter)
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

static std::expected<std::vector<Description::FullRange>, std::string>
handleFullRange (refimpls::FunctionType functype, const std::string &name)
{
  if (name == "normal")
    {
      switch (functype)
	{
	case refimpls::FunctionType::f32_f:
	case refimpls::FunctionType::f32_f_fp_fp:
	  return std::vector<Description::FullRange>{
	    Description::FullRange{ "positive normal (float)",
				    floatrange::Limits<float>::PlusNormalMin,
				    floatrange::Limits<float>::PlusNormalMax },
	    Description::FullRange{ "negative normal (float)",
				    floatrange::Limits<float>::NegNormalMin,
				    floatrange::Limits<float>::NegNormalMax },
	  };
	case refimpls::FunctionType::f64_f:
	case refimpls::FunctionType::f64_f_fp_fp:
	  return std::vector<Description::FullRange>{
	    Description::FullRange{
		"positive normal (double)",
		floatrange::Limits<double>::PlusNormalMin,
		floatrange::Limits<double>::PlusNormalMax },
	    Description::FullRange{ "negative normal (double)",
				    floatrange::Limits<double>::NegNormalMin,
				    floatrange::Limits<double>::NegNormalMax },
	  };
	default:
	  return std::unexpected (std::string ("invalid function type"));
	}
    }
  else if (name == "subnormal")
    {
      switch (functype)
	{
	case refimpls::FunctionType::f32_f:
	case refimpls::FunctionType::f32_f_fp_fp:
	  return std::vector<Description::FullRange>{
	    Description::FullRange{
		"positive subnormal (float)",
		floatrange::Limits<float>::PlusSubnormalMin,
		floatrange::Limits<float>::PlusSubnormalMax },
	    Description::FullRange{
		"negative subnormal (float)",
		floatrange::Limits<float>::NegSubnormalMin,
		floatrange::Limits<float>::NegSubnormalMax },

	  };
	case refimpls::FunctionType::f64_f:
	case refimpls::FunctionType::f64_f_fp_fp:
	  return std::vector<Description::FullRange>{
	    Description::FullRange{
		"positive subnormal (float)",
		floatrange::Limits<double>::PlusSubnormalMin,
		floatrange::Limits<double>::PlusSubnormalMax },
	    Description::FullRange{
		"negative subnormal (float)",
		floatrange::Limits<double>::NegSubnormalMin,
		floatrange::Limits<double>::NegSubnormalMax },

	  };
	default:
	  return std::unexpected (std::string ("invalid function type"));
	}
    }
  return std::unexpected (std::format ("invalid full range: {}", name));
}

template <typename F>
inline std::expected<F, std::string>
parseRange (const std::string &str)
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
  return floatrange::fromStr<F> (std::string (trimmed));
}

template <>
inline std::expected<long long int, std::string>
parseRange (const std::string &str)
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

static std::expected<Description::SampleType, std::string>
handle1Arg (refimpls::FunctionType functype, const std::string &start,
	    const std::string &end, uint64_t count)
{
  switch (functype)
    {
    case refimpls::FunctionType::f32_f:
      return Description::SampleType (Description::Sample1Arg<float>{
	  TRY (parseRange<float> (start)), TRY (parseRange<float> (end)),
	  count });
    case refimpls::FunctionType::f64_f:
      return Description::SampleType (Description::Sample1Arg<double>{
	  TRY (parseRange<double> (start)), TRY (parseRange<double> (end)),
	  count });
    default:
      std::unreachable ();
    }
}

static std::expected<Description::SampleType, std::string>
handle2Arg (refimpls::FunctionType functype, const std::string &start_x,
	    const std::string &end_x, const std::string &start_y,
	    const std::string &end_y, uint64_t count)
{
  switch (functype)
    {
    case refimpls::FunctionType::f32_f_f:
      return Description::SampleType (Description::Sample2Arg<float>{
	  TRY (parseRange<float> (start_x)), TRY (parseRange<float> (end_x)),
	  TRY (parseRange<float> (start_y)), TRY (parseRange<float> (end_y)),
	  count });
    case refimpls::FunctionType::f64_f_f:
      return Description::SampleType (Description::Sample2Arg<double>{
	  TRY (parseRange<double> (start_x)), TRY (parseRange<double> (end_x)),
	  TRY (parseRange<double> (start_y)), TRY (parseRange<double> (end_y)),
	  count });
    case refimpls::FunctionType::f32_f_lli:
      return Description::SampleType (Description::Sample2ArgLli<float>{
	  TRY (parseRange<float> (start_x)), TRY (parseRange<float> (end_x)),
	  TRY (parseRange<long long int> (start_y)),
	  TRY (parseRange<long long int> (end_y)), count });
    case refimpls::FunctionType::f64_f_lli:
      return Description::SampleType (Description::Sample2ArgLli<double>{
	  TRY (parseRange<double> (start_x)), TRY (parseRange<double> (end_x)),
	  TRY (parseRange<long long int> (start_y)),
	  TRY (parseRange<long long int> (end_y)), count });
    default:
      std::unreachable ();
    }
}

template <class... Ts> struct overloaded : Ts...
{
  using Ts::operator()...;
};

std::expected<void, std::string>
Description::parse (const std::string &fname)
{
  std::ifstream file (fname);

  nlohmann::json data;
  try
    {
      data = nlohmann::json::parse (file);
    }
  catch (nlohmann::json::parse_error &ex)
    {
      return std::unexpected (std::format (
	  "error parsing file {} at position {}", fname, ex.byte));
    }

  this->FunctionName = data["function"];
  auto functype = refimpls::getFunctionType (FunctionName);
  if (!functype)
    return std::unexpected (
	std::format ("invalid FunctionName: {}", FunctionName));

  if (data.contains ("full"))
    {
      auto fulldata
	  = SplitWithRanges (data["full"].get<std::string_view> (), ",");
      for (const auto &f : fulldata)
	{
	  auto full = handleFullRange (functype.value (), f);
	  if (!full)
	    return std::unexpected (full.error ());
	  auto fullv = full.value ();
	  this->Samples.insert (this->Samples.end (), fullv.begin (),
				fullv.end ());
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

	      auto sample = TRY (handle2Arg (
		  functype.value (), r["x"][0].template get<std::string> (),
		  r["x"][1].template get<std::string> (),
		  r["y"][0].template get<std::string> (),
		  r["y"][1].template get<std::string> (),
		  r["count"].get<uint64_t> ()));
	      this->Samples.push_back (sample);
	    }
	  else if (r.contains ("x"))
	    {
	      if (r["x"].size () != 2)
		return std::unexpected (std::format (
		    "invalid sample size: {} (expected 2)", r["x"].size ()));

	      auto sample = TRY (handle1Arg (
		  functype.value (), r["x"][0].template get<std::string> (),
		  r["x"][1].template get<std::string> (),
		  r["count"].get<uint64_t> ()));

	      this->Samples.push_back (sample);
	    }
	  else
	    return std::unexpected (
		std::format ("invalid sample definition {}", r.dump ()));
	}
    }
  else
    return std::unexpected (std::string ("no samples found"));

  return {};
}
