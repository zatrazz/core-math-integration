#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <ctype.h>
#include <stdbool.h>
#include <errno.h>
#include <float.h>

static inline bool
startswith (const char *str, const char *pre)
{
  while (isspace (*str))
    str++;

  size_t lenpre = strlen (pre);
  size_t lenstr = strlen (str);
  return lenstr < lenpre ? false : memcmp (pre, str, lenpre) == 0;
}

int main (int argc, char *argv[])
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

  return 0;
}
