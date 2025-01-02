#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <ctype.h>
#include <stdbool.h>
#include <errno.h>
#include <float.h>
#include <unistd.h>

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

static void
check_float (void)
{

  float max = FLT_MIN;
  float min = FLT_MAX;

  char *line = NULL;
  size_t linesz = 0;
  while (getline (&line, &linesz, stdin) != -1)
    {
      if (startswith (line, "#"))
	continue;

      char *endptr;
      float f = strtof (line, &endptr);

      if (f < min)
	min = f;
      if (f > max)
	max = f;
    }

  printf ("min=%a (%f %g)  max=%a (%f %g)\n",
	  min, min, min,
	  max, max, max);
}

static void
check_double (void)
{
  double max = DBL_MIN;
  double min = DBL_MAX;

  char *line = NULL;
  size_t linesz = 0;
  while (getline (&line, &linesz, stdin) != -1)
    {
      if (startswith (line, "#"))
	continue;

      char *endptr;
      double f = strtof (line, &endptr);

      if (f < min)
	min = f;
      if (f > max)
	max = f;
    }

  printf ("min=%la (%lf %lg)  max=%la (%lf %lg)\n",
	  min, min, min,
	  max, max, max);
}

int main (int argc, char *argv[])
{
  enum ftype { FLOAT, DOUBLE } type = FLOAT;

  int opt;
  while ((opt = getopt (argc, argv, "t:")) != -1)
    {
      switch (opt)
	{
	case 't':
	  if (strequal (optarg, "float"))
	    type = FLOAT;
	  else if (strequal (optarg, "double"))
	    type = DOUBLE;
	  else
	    {
	      fprintf (stderr, "error: invalid type (%s)\n", optarg);
	      exit (EXIT_FAILURE);
	    }
	  break;

	default:
	  fprintf (stderr, "usage: %s [-t type] input\n", argv[0]);
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
    case FLOAT:  check_float (); break;
    case DOUBLE: check_double (); break;
    }

  return 0;
}
