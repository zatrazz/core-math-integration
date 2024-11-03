#include <stdio.h>
#include <math.h>

float cr_lgammaf (float);

int main ()
{
#if 0
  float start = 0x1.0a6d2p+123;

  while (1)
    {
      start = nexttowardf (start, 0.0f);
      float r = cr_lgammaf (start);
      //printf ("start=%a (%a)\n", start, r);
      if (isfinite (r))
	break;
    }

  printf ("start=%a %f %g\n",
	  start, start, start);

  start = nexttowardf (start, INFINITY);
  float r = cr_lgammaf (start);
  printf ("nextt=%a %f %g | %a %f %g\n",
	  start, start, start,
	  r, r, r);
#else

  float x = cr_lgammaf (-0x1.78p+8);
  printf ("r=%a\n", x);
#endif
}
