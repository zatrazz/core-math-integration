#include <math.h>
#include <stdint.h>

// NB: this only works for finite numbers, which is suffice for
// the CORE-MATH usage.
float
roundevenf (float x)
{
  float y = roundf (x);
  if (fabs (x - y) == 0.5f)
    {
      union
      {
        float f;
        uint32_t i;
      } u = { y };
      union
      {
        float f;
        uint32_t i;
      } v = { y - copysignf (1.0f, x) };
      if (__builtin_ctzl (v.i) > __builtin_ctzl (u.i))
        y = v.f;
    }
  return y;
}
