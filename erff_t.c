#include <stdio.h>
#include <math.h>

extern float cr_erff (float);

int main ()
{
  float n = 0x5.ebed23725989cp+0;

  int i = 0;
  while (1)
    {
      float r = cr_erff (n);
      if (r < 1.0f)
	break;
      n = nexttowardf (n, 0.0f);
      i += 1;
    }
  printf ("n=%a | %f | %g (%d)\n", n, n, n, i);

  return 0;
}
