#ifndef _CXXCOMPAT_H

#include "config.h"

#ifdef HAVE_PRINT_HEADER
# include <print>
#else
#include <iostream>
#include <format>

namespace std {
template <typename... Args>
inline void
println (const std::format_string<Args...> fmt, Args &&...args)
{
  std::cout << std::vformat (fmt.get (), std::make_format_args (args...))
            << '\n';
}
};
#endif

#endif
