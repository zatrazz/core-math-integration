#ifndef REFIMPLS_H
#define REFIMPLS_H

#include <expected>
#include <functional>
#include <string_view>

namespace refimpls
{

typedef std::function<float (float)> univariate_binary32_t;
typedef std::function<float (float, int rnd)> univariate_binary32_ref_t;

typedef std::function<double (double)> univariate_binary64_t;
typedef std::function<double (double, int rnd)> univariate_binary64_ref_t;

typedef std::function<float (float, float)> bivariate_binary32_t;
typedef std::function<float (float, float, int rnd)>
    bivariate_binary32_ref_t;

typedef std::function<double (double, double)> bivariate_binary64_t;
typedef std::function<double (double, double, int rnd)>
    bivariate_binary64_ref_t;

enum class errors_t
{
  invalid_func
};

enum class func_type_t
{
  binary32_univariate,
  binary32_bivariate,
  binary64_univariate,
  binary64_bivariate,
};

std::expected<func_type_t, errors_t> get_func_type (const std::string_view &);

std::expected<std::pair<univariate_binary32_t, univariate_binary32_ref_t>,
              errors_t>
get_univariate_binary32 (const std::string_view &, bool coremath = false);

std::expected<std::pair<bivariate_binary32_t, bivariate_binary32_ref_t>,
              errors_t>
get_bivariate_binary32 (const std::string_view &, bool coremath = false);

std::expected<std::pair<univariate_binary64_t, univariate_binary64_ref_t>,
              errors_t>
get_univariate_binary64 (const std::string_view &, bool coremath = false);

std::expected<std::pair<bivariate_binary64_t, bivariate_binary64_ref_t>,
              errors_t>
get_bivariate_binary64 (const std::string_view &, bool coremath = false);

} // refimpls

#endif
