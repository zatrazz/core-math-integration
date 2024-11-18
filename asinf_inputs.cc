#include <iostream>
#include <iomanip>
#include <string>
#include <map>
#include <random>
#include <charconv>
#include <string_view>
#include <cassert>

extern "C" {
    float cr_asinf (float);
};

typedef std::numeric_limits<float> binary32;

int
main (int argc, char *argv[])
{
  constexpr const int n = 1000;

  std::random_device rd;

  std::mt19937 e2 (rd ());
  //std::knuth_b e2(rd());
  //std::default_random_engine e2(rd()) ;

  const struct
  {
    const char *desc;
    const char *title;
    const float start;
    const float end;
    const int n;
  } ranges[] = {
    {
      "Random inputs in [-1,1]",
      "random",
      -0x1p+0,
       0x1p+0,
      2705,
    },
  };


  std::printf ("## args: float\n");
  std::printf ("## ret: float\n");
  std::printf ("## includes: math.h\n");
  for (auto range : ranges)
    {
      std::uniform_real_distribution<float> dist(range.start, range.end);

#if 0
      float max = dist1.max();
      float min = dist1.min();
      std::printf ("%a (%f) %a (%f)\n", min, min, max, max);
#endif

      std::printf ("# %s\n", range.desc);
      std::printf ("## name: workload-%s\n", range.title);
      for (int i = 0; i < range.n; i++)
	{
	  float f = dist (e2);
	  float r = cr_asinf (f);
	  assert (std::isfinite (r));
	  std::printf ("%a\n", f);
	}
    }
}
