#include <cassert>
#include <random>

extern "C" {
    float cr_cbrtf (float);
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
    const float start;
    const float end;
    const int n;
  } ranges[] = {
    {
      "Random inputs in [1,8]",
      "random-1-8",
      0x1p+0,
      0x1p+3,
      1000,
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
	  float r = cr_cbrtf (f);
	  assert (std::isfinite (r));
	  std::printf ("%a\n", f);
	}
    }
}
