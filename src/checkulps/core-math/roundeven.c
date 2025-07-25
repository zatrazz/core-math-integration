#include <math.h>
#include <stdint.h>

// NB: this only works for finite numbers, which is suffice for
// the CORE-MATH usage.
double
roundeven (double x)
{
  double y = round (x);
  if (fabs (x - y) == 0.5)
    {
      union
      {
        double f;
        uint64_t i;
      } u = { y };
      union
      {
        double f;
        uint64_t i;
      } v = { y - copysign (1.0, x) };
      if (__builtin_ctzll (v.i) > __builtin_ctzll (u.i))
        y = v.f;
    }
  return y;
}
