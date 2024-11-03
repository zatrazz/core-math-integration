#include <cstdio>
#include <random>
#include <ieee754.h>
#include <cassert>

typedef std::numeric_limits<float> binary32;

int
main (int argc, char *argv[])
{
  constexpr const int n = 1000;

  std::random_device rd;
  std::default_random_engine e (rd()) ;

  std::printf ("## args: float\n");
  std::printf ("## ret: float\n");
  std::printf ("## include: math.h\n");
  std::printf ("# Random inputs x*2^e where x is random in [1/2,1] and e in [-29,127]\n");
  std::printf ("## name: workload-random\n");

  std::uniform_real_distribution<float> dist_x (0x1p-1, 0x1p+0);
  std::uniform_int_distribution<int> dist_e (-29, 127);

  for (auto i = 0; i < n; i++)
    {
#if 0
      union ieee754_float f =
      {
	.ieee =
	{
	  .negative = 0;  // as per log1p-inputs
	  .exponent = dist_e (e),
	  .mantissa = dist_m (A

	}
      };
#endif
      float x = dist_x (e);
      float exp = powf (2.0, (float) dist_e (e));
      float r = x * exp;

      assert (std::isfinite (std::log1pf (r)));

      std::printf ("%a\n", r);
    }
}
