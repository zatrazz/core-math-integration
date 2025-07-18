#ifndef REFIMPLS_H
#define REFIMPLS_H

#include <expected>
#include <functional>
#include <string_view>

namespace refimpls
{

template <typename F> using univariate_t = std::function<F (F)>;
template <typename F> using univariate_ref_t = std::function<F (F, int)>;

template <typename F> using bivariate_t = std::function<F (F, F)>;
template <typename F> using bivariate_ref_t = std::function<F (F, F, int)>;

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

template <class F> void init_ref_func ();

std::expected<func_type_t, errors_t> get_func_type (const std::string_view &);

template <typename F>
std::expected<std::pair<univariate_t<F>, univariate_ref_t<F> >, errors_t>
get_univariate (const std::string_view &, bool coremath = false);

template <typename F>
std::expected<std::pair<bivariate_t<F>, bivariate_ref_t<F> >, errors_t>
get_bivariate (const std::string_view &, bool coremath = false);

} // refimpls

#endif
