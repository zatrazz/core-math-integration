//
// Copyright (c) Adhemerval Zanella. All rights reserved.
// Licensed under the MIT license. See LICENSE file in the project root for
// details.
//

#ifndef _IOHELPER_H
#define _IOHELPER_H

#include <iostream>
#include <chrono>
#include "cxxcompat.h"

namespace iohelper {

//
// Helper output functions.
//

static inline void
println_ts (const std::string &str)
{
  auto now = std::chrono::system_clock::now ();
  auto seconds = std::chrono::floor<std::chrono::seconds> (now);
  std::println ("[{:%Y-%m-%d %H:%M:%S}] {}", seconds, str);
}

template <typename... Args>
static inline void
println_ts (std::format_string<Args...> fmt, Args &&...args)
{
  auto now = std::chrono::system_clock::now ();
  auto seconds = std::chrono::floor<std::chrono::seconds> (now);
  std::print ("[{:%Y-%m-%d %H:%M:%S}] ", seconds);
  std::println (fmt, std::forward<Args> (args)...);
}

template <typename... Args>
[[noreturn]] inline void
error (const std::format_string<Args...> fmt, Args &&...args)
{
  std::cerr << "error: "
            << std::vformat (fmt.get (), std::make_format_args (args...))
            << '\n';
  std::exit (EXIT_FAILURE);
}

[[noreturn]] static inline void
error (const std::string &str)
{
  std::cerr << "error: " << str << '\n';
  std::exit (EXIT_FAILURE);
}

}

#endif
