#include <iostream>
#include <iomanip>
#include <string>
#include <map>
#include <random>
#include <charconv>
#include <string_view>
#include <cassert>

extern "C" {
    float cr_acoshf (float);
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
      "Random inputs in the range [1,21]",
      "random-1-21",
      0x1p+0,
      0x1.5p+4,
      1000,
    },
#if 0
    {
      "first polinomial used by CORE-MATH [1.0, 1.202000)",
      "core-math1",
      0x1p+0,
      0x1.33b646p+0,
    },
    {
      "second polinomial used by CORE-MATH (1.202000, FLT_MAX)",
      "core-math2",
      0x1.33b648p+0,
      std::numeric_limits<float>::max(),
    },
    {
      "approximation with log1pf [1.0, 2.0]",
      "approximation1",
      0x1p+0,
      0x1p+1,
    },
    {
      "approximation with logf + sqrt (2.0, 2**28]",
      "approximation2",
      0x1.000002p+1,
      0x1p+28,
    },
    {
      "approximation with logf (2**28)",
      "approximation3",
      0x1.000002p+28,
      std::numeric_limits<float>::max(),
    },
#endif
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
	  float r = cr_acoshf (f);
	  assert (std::isfinite (r));
	  std::printf ("%a\n", f);
	}
    }
}
