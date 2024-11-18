#include <cstdio>
#include <random>
#include <ieee754.h>
#include <cassert>
#include <numbers>

typedef std::numeric_limits<float> binary32;

extern "C" {
    float cr_tanf (float);
};

int
main (int argc, char *argv[])
{
  constexpr const int n = 3000;

  std::random_device rd;
  std::default_random_engine e (rd()) ;

  std::printf ("## args: float\n");
  std::printf ("## ret: float\n");
  std::printf ("## includes: math.h\n");

  {
    std::printf ("# Random inputs in [-pi, pi]\n");
    std::printf ("## name: workload-pi-pi\n");
  
    constexpr const float range_start = -std::numbers::pi_v<float>;
    constexpr const float range_end =    std::numbers::pi_v<float>;
  
    std::uniform_real_distribution<float> dist_x (range_start, range_end);
  
    for (auto i = 0; i < n; i++)
      {
        float x = dist_x (e);
        float r = cr_tanf (x);
  
        assert (std::isfinite (r));
  
        std::printf ("%a\n", x);
      }
  }

  {
    std::printf ("# Random inputs with large exponent\n");
    std::printf ("## name: workload-rbig\n");

    constexpr const float x_start = 0x1p-1;
    constexpr const float x_end   = 0x1p+0;
    constexpr const int   e_start = 28;
    constexpr const int   e_end   = (255 - 1) - 127;

    std::uniform_real_distribution<float> dist_x (x_start, x_end);
    std::uniform_int_distribution<int> dist_e (e_start, e_end);

    for (auto i = 0; i < n; i++)
      {
	float x_rand = dist_x (e);
	int   e_rand = dist_e (e);

	float r_rand = ldexpf (x_rand, e_rand);

	float r_tanf = cr_tanf (r_rand);
	//assert (std::isfinite (r_tanf));

	std::printf ("%a\n", r_rand);
      }
  }
}
