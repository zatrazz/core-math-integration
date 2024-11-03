#include <cstdio>
#include <random>
#include <cassert>

extern "C"
{
  float cr_powf (float, float);
  float cr_log10f (float);
};

typedef std::numeric_limits<float> binary32;

int
main (int argc, char *argv[])
{
  constexpr const int n = 1000;
  constexpr const float x_s = 0x1p-1;
  constexpr const float x_e = 0x1p+0;
  constexpr const int e_s = -126;
  constexpr const int e_e = 127;;

  std::random_device rd;
  std::default_random_engine e (rd()) ;

  std::printf ("## args: float\n");
  std::printf ("## ret: float\n");
  std::printf ("## includes: math.h\n");
  std::printf ("# Random inputs x*2^e where x is random "
	       "in [%a,%a] and e in [%d,%d]\n",
    x_s, x_e, e_s, e_e);
  std::printf ("## name: workload-random\n");

  std::uniform_real_distribution<float> dist_x (x_s, x_e);
  std::uniform_int_distribution<int> dist_e (e_s, e_e);

  for (auto i = 0; i < n; i++)
    {
      float x = dist_x (e);
      float exp = cr_powf (2.0f, (float) dist_e (e));
      float r = x * exp;

      assert (std::isfinite (cr_log10f (r)));

      std::printf ("%a\n", r);
    }
}
