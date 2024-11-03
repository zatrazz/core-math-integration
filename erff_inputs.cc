#include <cstdio>
#include <random>
#include <cassert>

extern "C"
{
  float cr_erff (float);
};

typedef std::numeric_limits<float> binary32;

int
main (int argc, char *argv[])
{
  constexpr const float x_start = 0x0p+0;
  constexpr const float x_end   = 0x1.f5a888p+1;
  constexpr const int n = 794 - 5;

  std::random_device rd;
  std::default_random_engine e (rd()) ;

  std::printf ("## args: float\n");
  std::printf ("## ret: float\n");
  std::printf ("## includes: math.h\n");
  std::printf ("# Random inputs in [%1.0f,b=%a]\n", x_start, x_end);
  std::printf ("# where b in the smallest number such that erff(b) ");
  std::printf ("rounds to 1 (to nearest)\n");
  std::printf ("## name: workload-random\n");

  std::uniform_real_distribution<float> dist_x (x_start, x_end);

  for (auto i = 0; i < n; i++)
    {
      float x = dist_x (e);
      float r_cr = cr_erff (x);

      assert (std::isfinite (r_cr));
      assert (r_cr < 1.0f);

      std::printf ("%a\n", x);
    }
}
