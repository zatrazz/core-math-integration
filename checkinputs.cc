#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <ctype.h>
#include <stdbool.h>
#include <errno.h>
#include <float.h>
#include <unistd.h>

#include <string>
#include <iostream>
#include <limits>
#include <charconv>
#include <format>
#include <system_error>

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

template<typename... Args>
inline void println(const std::format_string<Args...> fmt, Args&&... args)
{
  std::cout << std::vformat(fmt.get(), std::make_format_args(args...)) << '\n';
}

template <typename ftype, ftype (strtotype)(const char *, char **)>
void check_type (void)
{
  ftype max = std::numeric_limits<ftype>::min();
  ftype min = std::numeric_limits<ftype>::max();

  std::string str;
  while(std::getline (std::cin, str))
    {
      if (str.rfind ('#', 0) == 0)
	continue;

      ftype f = strtotype (str.c_str(), NULL);
      if (f < min)
	min = f;
      if (f > max)
	max = f;
    }

  println ("min={:a} ({:g})  max={:a} {:g})",
	   min, min,
	   max, max);
}

template <typename ftype, ftype (strtotype)(const char *, char **)>
void check_type_bivariate (void)
{
  ftype max = std::numeric_limits<ftype>::min();
  ftype min = std::numeric_limits<ftype>::max();

  std::string str;
  while(std::getline (std::cin, str))
    {
      if (str.rfind ('#', 0) == 0)
	continue;

      char *end;
      ftype f0 = strtotype (str.c_str(), &end);
      if (*end != ',')
	continue;

      ftype f1 = strtotype (end + 1, NULL);
      
      if (f0 < min)
	min = f0;
      if (f0 > max)
	max = f0;

      if (f1 < min)
	min = f1;
      if (f1 > max)
	max = f1;
    }

  println ("min={:a} ({:g})  max={:a} ({:g})",
	   min, min,
	   max, max);
}

int main (int argc, char *argv[])
{
  enum ftype { FLOAT, DOUBLE, LDOUBLE } type = FLOAT;
  bool bivariate = false;

  int opt;
  while ((opt = getopt (argc, argv, "t:b")) != -1)
    {
      switch (opt)
	{
	case 't':
	  if (strequal (optarg, "binary32"))
	    type = FLOAT;
	  else if (strequal (optarg, "binary64"))
	    type = DOUBLE;
	  else if (strequal (optarg, "binary80"))
	    type = LDOUBLE;
	  else
	    {
	      fprintf (stderr, "error: invalid type (%s)\n", optarg);
	      exit (EXIT_FAILURE);
	    }
	  break;
	
	case 'b':
	  bivariate = true;
	  break;

	default:
	  fprintf (stderr, "usage: %s [-t type] [-b] input\n", argv[0]);
	  exit (EXIT_FAILURE);
	}
    }

#if 0
  if (optind >= argc)
    {
      fprintf (stderr, "error: expected argument after options\n");
      exit (EXIT_FAILURE);
    }
#endif

  switch (type)
    {
    case FLOAT:  
      if (bivariate)
	check_type_bivariate<float, std::strtof> ();
      else
	check_type<float, std::strtof> ();
      break;

    case DOUBLE:
      if (bivariate)
	check_type_bivariate<double, std::strtod> ();
      else
	check_type<double, std::strtod> ();
      break;

    case LDOUBLE:
      if (bivariate)
	check_type_bivariate<long double, std::strtold> ();
      else
	check_type<long double, std::strtold> ();
      break;
    }

  return 0;
}
