#include <algorithm>
#include <cctype>
#include <cstring>
#include <format>
#include <iostream>
#include <limits>
#include <locale>
#include <string>

#include <boost/program_options.hpp>

template <typename... Args>
inline void
println (const std::format_string<Args...> fmt, Args &&...args)
{
  std::cout << std::vformat (fmt.get (), std::make_format_args (args...))
            << '\n';
}

template <typename T> struct ftypeinfo
{
  static std::string name ();
  static T fromstring (const std::string &);
};

template <> struct ftypeinfo<float>
{
  using type = float;
  static std::string
  name ()
  {
    return "float";
  }
  static float
  fromstring (const std::string &s)
  {
    return std::stof (s);
  }
};

template <> struct ftypeinfo<double>
{
  using type = double;
  static std::string
  name ()
  {
    return "double";
  }
  static double
  fromstring (const std::string &s)
  {
    return std::stod (s);
  }
};

template <> struct ftypeinfo<long double>
{
  using type = long double;
  static std::string
  name ()
  {
    return "long double";
  }
  static long double
  fromstring (const std::string &s)
  {
    return std::stold (s);
  }
};

static inline bool
startswith (const char *str, const char *pre)
{
  while (isspace (*str))
    str++;

  size_t lenpre = strlen (pre);
  size_t lenstr = strlen (str);
  return lenstr < lenpre ? false : memcmp (pre, str, lenpre) == 0;
}

static inline bool
strequal (const char *str1, const char *str2)
{
  return strcmp (str1, str2) == 0;
}

// trim from start (in place)
inline void
ltrim (std::string &s)
{
  s.erase (s.begin (),
           std::find_if (s.begin (), s.end (), [] (unsigned char ch) {
             return !std::isspace (ch);
           }));
}

// trim from end (in place)
inline void
rtrim (std::string &s)
{
  s.erase (std::find_if (s.rbegin (), s.rend (),
                         [] (unsigned char ch) { return !std::isspace (ch); })
               .base (),
           s.end ());
}

inline void
trim (std::string &s)
{
  rtrim (s);
  ltrim (s);
}

template <typename ftypeinfo>
static void
check_type (void)
{
  using ftype = typename ftypeinfo::type;

  std::vector<ftype> numbers;
  std::string str;
  while (std::getline (std::cin, str))
    {
      if (str.rfind ('#', 0) == 0)
        continue;

      trim (str);
      if (str.empty ())
        continue;

      numbers.push_back (ftypeinfo::fromstring (str));
    }
  if (numbers.empty ())
    return;

  auto minmaxt = std::minmax_element (numbers.begin (), numbers.end ());

  println ("min={0:a} ({0:g})  max={1:a} ({1:g})", *minmaxt.first,
           *minmaxt.second);
}

template <typename ftype, ftype (strtotype) (const char *, char **)>
void
check_type_bivariate (void)
{
  std::vector<ftype> numbers0;
  std::vector<ftype> numbers1;
  std::string str;
  while (std::getline (std::cin, str))
    {
      if (str.rfind ('#', 0) == 0)
        continue;

      trim (str);
      if (str.empty ())
        continue;

      char *end;
      ftype f0 = strtotype (str.c_str (), &end);
      if (*end != ',')
        continue;
      numbers0.push_back (f0);

      numbers1.push_back (strtotype (end + 1, NULL));
    }
  if (numbers0.empty () || numbers1.empty ())
    return;

  auto minmaxt0 = std::minmax_element (numbers0.begin (), numbers0.end ());
  auto minmaxt1 = std::minmax_element (numbers1.begin (), numbers1.end ());

  println ("min={0:a}, {1:a} ({0:g}, {1:g})  max={2:a}, {3:a} ({2:g}, {3:g})",
           *minmaxt0.first, *minmaxt0.second, *minmaxt1.first,
           *minmaxt1.second);
}

int
main (int argc, char *argv[])
{
  namespace po = boost::program_options;
  po::options_description desc{ "options" };
  desc.add_options () ("help,h", "help") (
      "type,t", po::value<std::string> ()->default_value ("binary32"),
      "type to use") ("bivariate,b",
                      po::bool_switch ()->default_value (false));

  po::variables_map vm;
  po::store (po::parse_command_line (argc, argv, desc), vm);
  po::notify (vm);

  if (vm.count ("help"))
    {
      std::cout << desc << "\n";
      return 1;
    }

  std::string type = vm["type"].as<std::string> ();
  bool bivariate = vm["bivariate"].as<bool> ();

  if (type == "binary32")
    {
      if (bivariate)
        check_type_bivariate<float, std::strtof> ();
      else
        check_type<ftypeinfo<float> > ();
    }
  else if (type == "binary64")
    {
      if (bivariate)
        check_type_bivariate<double, std::strtod> ();
      else
        check_type<ftypeinfo<double> > ();
    }
  else if (type == "binary80")
    {
      if (bivariate)
        check_type_bivariate<long double, std::strtold> ();
      else
        check_type<ftypeinfo<long double> > ();
    }

  return 0;
}
