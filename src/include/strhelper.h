#ifndef _STRHELPER_H
#define _STRHELPER_H

namespace strhelper
{

static inline std::string_view
trim (const std::string &str)
{
  std::string_view sv{ str };

  size_t start = sv.find_first_not_of (" \t\n\r\f\v");
  if (start == std::string_view::npos)
    {
      return std::string_view{};
    }

  size_t end = sv.find_last_not_of (" \t\n\r\f\v");
  return sv.substr (start, end - start + 1);
}

} // namespace strhelper

#endif
