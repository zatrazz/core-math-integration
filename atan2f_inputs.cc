#include <iostream>
#include <iomanip>
#include <string>
#include <map>
#include <random>
#include <charconv>
#include <string_view>
#include <cassert>

extern "C" {
    float cr_atan2f (float, float);
};

typedef std::numeric_limits<float> binary32;

int
main (int argc, char *argv[])
{
  constexpr const int n = 2000;

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
  } ranges[] = {
    {
      "Random x,y inputs in the range [-10,10]",
      "random",
      -0x1.4p+3,
       0x1.4p+3,
    },
  };


  std::printf ("## args: float:float\n");
  std::printf ("## ret: float\n");
  std::printf ("## includes: math.h\n");
  for (auto range : ranges)
    {
      std::uniform_real_distribution<float> dist(range.start, range.end);

      std::printf ("# %s\n", range.desc);
      std::printf ("## name: workload-%s\n", range.title);
      for (int i = 0; i < n; i++)
	{
	  float x = dist (e2);
	  float y = dist (e2);
	  float r = cr_atan2f (x, y);
	  assert (std::isfinite (r));
	  std::printf ("%a, %a\n", x, y);
	}
    }
}
