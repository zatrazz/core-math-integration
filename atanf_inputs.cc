#include <iostream>
#include <iomanip>
#include <string>
#include <map>
#include <random>
#include <charconv>
#include <string_view>
#include <cassert>

extern "C" {
    float cr_atanf (float);
};

typedef std::numeric_limits<float> binary32;

int
main (int argc, char *argv[])
{
#if 0
  if (argc < 2)
    {
      std::printf ("usage: %s <n>\n", argv[0]);
      std::exit (EXIT_SUCCESS);
    }

  int n;
  {
    std::string_view str(argv[1]);
    auto [ptr, ec] = std::from_chars (str.data(), str.data() + str.size(), n);
    if (ec != std::errc())
      {
        std::printf ("error: invalid number (%s)\n", argv[1]);
        std::exit (EXIT_FAILURE);
      }
  }
#endif

  constexpr const int n = 2000;

  std::random_device rd;

  //std::mt19937 e2 (rd ());
  //std::knuth_b e2(rd());
  std::default_random_engine e2(rd()) ;

  const struct
  {
    const char *desc;
    const char *title;
    float start;
    float end;
  } ranges[] = {
    {
      "Random inputs in the range [-10,10]",
      "core-math1",
      -0x1.4p+3,
       0x1.4p+3,
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
      for (int i = 0; i < n; i++)
	{
	  float f = dist (e2);
	  float r = cr_atanf (f);
	  assert (std::isfinite (r));
	  assert (f >= range.start);
	  assert (f <= range.end);
	  std::printf ("%a\n", f);
	}
    }
}
