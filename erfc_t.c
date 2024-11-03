#include <stdio.h>
#include <math.h>
#include <stdbool.h>

extern float cr_erfcf (float);

int main ()
{
#if 1
  float n = 1.0;
  float r;

  int i = 0;
  //bool print = false;

  while (1)
    {
      float next = nexttowardf (n, INFINITY);
      //if (print)
	//printf ("testing=(%a,%f,%g): ", next, next, next);
      r = cr_erfcf (next);
      if (r == 0.0f)
	break;
      //if (print)
        //printf ("(%a,%f,%g) (%d)\n", r, r, r, i);
      n = next;
      i += 1;
      //if (i == (27319803 - 1024))
	//print = false;
    }

  r = cr_erfcf (n);
  printf ("(%a,%f,%g)=(%a,%f,%g) (%d)\n",
	  n, n, n,
	  r, r, r,
	  i);

  n = nexttowardf (n, INFINITY);
  r = cr_erfcf (n);
  printf ("(%a,%f,%g)=(%a,%f,%g) (%d)\n",
	  n, n, n,
	  r, r, r,
	  i);


#else
  printf ("%a -> %a\n",
	  0x1.41bbf6p+3f,
	  cr_erfcf (0x1.41bbf6p+3f));
  printf ("%a -> %a\n",
	  nexttowardf (0x1.41bbf6p+3f, INFINITY),
	  cr_erfcf (nexttowardf (0x1.41bbf6p+3f, INFINITY)));
#endif


  return 0;
}
