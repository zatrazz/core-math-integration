#include <cassert>
#include <charconv>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <random>
#include <string>

int
main (int argc, char *argv[])
{
  if (argc < 2)
    {
      std::printf ("usage: %s <n>\n", argv[0]);
      exit (EXIT_SUCCESS);
    }

  int n;
  {
    std::string_view str = argv[1];
    auto [ptr, ec] = std::from_chars(str.data(), str.data() + str.size(), n);
    if (ec != std::errc())
      {
	std::fprintf (stderr, "error: invalid number\n");
	exit (EXIT_FAILURE);
      }
  }

  std::random_device rd;

  std::mt19937 gen (rd ());

  constexpr float start = -0x1.cfdadap+3;
  constexpr float end   =  0x1.62e42ep+6;

  std::uniform_real_distribution <> dist (start, end);
  //std::normal_distribution<> dist(start, end);
  //std::student_t_distribution<> dist(end);
  //std::poisson_distribution<> dist(start);
  //std::extreme_value_distribution<> dist(start, end);

  std::printf ("## args: float\n");
  std::printf ("## ret: float\n");
  std::printf ("## includes: math.h\n");
  std::printf ("# Random inputs in [a=%f,b=%f]\n", start, end);
  std::printf ("# where a is the smallest number such that expm1f does not round to -1\n");
  std::printf ("# and b is the smallest number such that expm1 rounds to +Inf (to nearest)\n");
  for (int i = 0; i < n; i++)
    {
      float n = dist (gen);
      float r = expm1f (n);
      assert (r != -1.0f);
      assert (!std::isinf (r));
      std::printf ("%a\n", n);
    }
}
