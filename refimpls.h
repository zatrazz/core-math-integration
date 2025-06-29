#ifndef REFIMPLS_H
#define REFIMPLS_H

#include <functional>
#include <expected>
#include <string_view>

namespace refimpls
{

typedef std::function<double(double)> univariate_t;

enum class errors_t
{
  invalid_func
};

std::expected<univariate_t, errors_t>
get_univariate (const std::string_view&, bool coremath = false);

std::expected<univariate_t, errors_t>
get_univariate_ref (const std::string_view&, int rnd);

} // refimpls

#endif
