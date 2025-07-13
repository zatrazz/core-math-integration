#ifndef REFIMPLS_H
#define REFIMPLS_H

#include <expected>
#include <functional>
#include <string_view>

namespace refimpls
{

typedef std::function<double (double)> univariate_t;
typedef std::function<double (double, int rnd)> univariate_ref_t;

typedef std::function<double (double, double)> bivariate_t;
typedef std::function<double (double, double, int rnd)> bivariate_ref_t;

enum class errors_t
{
  invalid_func
};

enum class func_type_t
{
  univariate,
  bivariate,
};

std::expected<func_type_t, errors_t> get_func_type (const std::string_view &);

std::expected<std::pair<univariate_t, univariate_ref_t>, errors_t>
get_univariate (const std::string_view &, bool coremath = false);

std::expected<std::pair<bivariate_t, bivariate_ref_t>, errors_t>
get_bivariate (const std::string_view &, bool coremath = false);

} // refimpls

#endif
