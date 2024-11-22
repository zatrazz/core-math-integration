#include <iostream>
#include <iomanip>
#include <string>
#include <map>
#include <random>
#include <charconv>
#include <string_view>
#include <cassert>

extern "C" {
    float cr_atanhf (float);
};

typedef std::numeric_limits<float> binary32;

int
main (int argc, char *argv[])
{
  std::random_device rd;

  std::mt19937 e2 (rd ());
  //std::knuth_b e2(rd());
  //std::default_random_engine e2(rd()) ;

  const struct
  {
    const char *desc;
    const char *title;
    float start;
    float end;
    int n;
  } ranges[] = {
    {
      "Random inputs in the range [-1,1]",
      "random",
      -1.0f,
       1.0f,
      2000,
    },
  };


  std::printf ("## args: float\n");
  std::printf ("## ret: float\n");
  std::printf ("## includes: math.h\n");
  for (auto range : ranges)
    {
      std::uniform_real_distribution<float> dist(range.start, range.end);

      std::printf ("# %s\n", range.desc);
      std::printf ("## name: workload-%s\n", range.title);
      for (int i = 0; i < range.n; i++)
	{
	  float f = dist (e2);
	  float r = cr_atanhf (f);
	  assert (std::isfinite (r));
	  std::printf ("%a\n", f);
	}
    }
}
