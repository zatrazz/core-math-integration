//
// Copyright (c) Adhemerval Zanella. All rights reserved.
// Licensed under the MIT license. See LICENSE file in the project root for
// details.
//

#ifndef _STRHELPER_H
#define _STRHELPER_H

#include <vector>
#include <ranges>

namespace strhelper
{

static inline std::string_view
trim (const std::string &str)
{
  std::string_view sv{ str };

  size_t start = sv.find_first_not_of (" \t\n\r\f\v");
  if (start == std::string_view::npos)
    return std::string_view{};

  size_t end = sv.find_last_not_of (" \t\n\r\f\v");
  return sv.substr (start, end - start + 1);
}

// trim from start (in place)
inline std::string &
ltrim (std::string &s)
{
  s.erase (s.begin (),
	   std::find_if (s.begin (), s.end (), [] (unsigned char ch) {
	     return !std::isspace (ch);
	   }));
  return s;
}

// trim from end (in place)
inline std::string &
rtrim (std::string &s)
{
  s.erase (std::find_if (s.rbegin (), s.rend (),
			 [] (unsigned char ch) { return !std::isspace (ch); })
	       .base (),
	   s.end ());
  return s;
}

inline std::string &
trim (std::string &s)
{
  rtrim (s);
  ltrim (s);

  return s;
}

static constexpr std::vector<std::string>
split_with_ranges (const std::string_view &s, std::string_view delimiter)
{
  std::vector<std::string> tokens;
  for (const auto &subrange : s | std::views::split (delimiter))
    tokens.push_back (std::string (subrange.begin (), subrange.end ()));
  return tokens;
}

} // namespace strhelper

#endif
