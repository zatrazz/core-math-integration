#include <iostream>
#include <iomanip>
#include <string>
#include <map>
#include <random>
#include <charconv>
#include <string_view>
#include <cassert>

extern "C" {
    float cr_asinhf (float);
};

typedef std::numeric_limits<float> binary32;

int
main (int argc, char *argv[])
{
  std::random_device rd;

  std::mt19937 e2 (rd ());
  //std::knuth_b e2(rd());
  //std::default_random_engine e2(rd()) ;

  std::printf ("## args: float\n");
  std::printf ("## ret: float\n");
  std::printf ("## includes: math.h\n");

  const struct
  {
    const char *desc;
    const char *title;
    float start;
    float end;
    int n;
  } ranges[] = {
    {
      "Random inputs in the range [-10,10]",
      "random",
      -0x1.4p+3,
       0x1.4p+3,
      2000,
    },
  };

  for (auto range : ranges)
    {
      std::uniform_real_distribution<float> dist(range.start, range.end);
      std::printf ("# %s\n", range.desc);
      std::printf ("## name: workload-%s\n", range.title);
      for (int i = 0; i < range.n; i++)
	{
	  float f = dist (e2);
	  float r = cr_asinhf (f);
	  assert (std::isfinite (r));
	  std::printf ("%a\n", f);
	}
    }

#if 0
  { 
    std::printf ("# random floats in [2^(e-1),2^e) "
		 "random exponent e in [-13,14]\n");
    std::printf ("## name: workload-random\n");

    std::uniform_int_distribution<int> e_dist(-13,14);
    for (int i = 0; i < 2000; i++)
      {
	int e = e_dist (e2);
	float start = std::ldexp (2.0f, e - 1);
	float end = std::ldexp (2.0f, e);

	std::uniform_real_distribution<float> dist_n (start, end);
	float n = dist_n (e2);
	assert (std::isfinite (n));
	std::printf ("%a\n", n);
      }
  }
#endif
}
