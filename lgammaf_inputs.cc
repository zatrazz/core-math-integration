#include <cstdio>
#include <random>
#include <cassert>

extern "C"
{
  float cr_powf (float, float);
  float cr_lgammaf (float);
};

typedef std::numeric_limits<float> binary32;

int
main (int argc, char *argv[])
{
  std::random_device rd;
  std::default_random_engine e (rd()) ;

  std::printf ("## args: float\n");
  std::printf ("## ret: float\n");
  std::printf ("## includes: math.h\n");

#if 0
  {
    constexpr const float x_start = 0x1p-1;
    constexpr const float x_end   = 0x1p+0;
    constexpr const int   e_start = -126;
    constexpr const int   e_end   = 127;
    constexpr const int   n       = 1000;

    std::printf ("# Random inputs x*2^e where x is random in [%1.1f,%1.0f] ",
		 x_start, x_end);
    std::printf ("and e in [%d,%d] (and result is a finite number)\n",
		 e_start, e_end);
    std::printf ("## name: workload-random\n");

    std::uniform_real_distribution<float> dist_x (x_start, x_end);
    std::uniform_int_distribution<int> dist_e (e_start, e_end);

    for (auto i = 0; i < n; i++)
      {
	float x_rand = dist_x (e);
	int   e_rand = dist_e (e);

	float r_rand = cr_powf (x_rand, (float) e_rand);

	float r_lgamma = cr_lgammaf (r_rand);
     
	if (!std::isfinite (r_lgamma))
	  continue;

	std::printf ("%a\n", r_rand);
      }
  }
#endif

  {
    constexpr const float x_start = -0x1.4p+4;
    constexpr const float x_end   = 0x1.4p+4;
    constexpr const int   n       = 1000;

    std::printf ("# Random inputs in the range [-20.0,20.0]\n");
    std::printf ("## name: workload-random-m20-p20\n");


    std::uniform_real_distribution<float> dist_x (x_start, x_end);
    for (auto i = 0; i < n; i++)
      {
	float x_rand = dist_x (e);

	float r_lgamma = cr_lgammaf (x_rand);
     
	assert (std::isfinite (r_lgamma));
	assert (x_rand <= 20);
	assert (-20 <= x_rand);

	std::printf ("%a\n", x_rand);
      }
  }
}
