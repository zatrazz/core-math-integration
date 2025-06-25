#include <algorithm>
#include <format>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <numbers>
#include <random>

#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include <climits>
#include <fenv.h>
#include <omp.h>

#include "refimpls.h"

using namespace refimpls;

template <typename RNG>
RNG seed_rng (size_t state_size = RNG::state_size)
{
  constexpr auto word_size = RNG::word_size / CHAR_BIT;
  constexpr auto rd_call_coefficient = word_size / sizeof (unsigned int);

  std::vector<unsigned int> data(state_size * rd_call_coefficient);
  std::generate(data.begin(), data.end(), std::random_device{});
  std::seed_seq seq(data.begin(), data.end());

  RNG mt(seq);
  return mt;
}

template<typename... Args>
inline void println(const std::format_string<Args...> fmt, Args&&... args)
{
  std::cout << std::vformat(fmt.get(), std::make_format_args(args...)) << '\n';
}

/* Returns the size of an ulp for VALUE.  */
template <typename F> F ulp (F value)
{
  F ulp;

  switch (std::fpclassify (value))
    {
      case FP_ZERO:
        /* We compute the distance to the next FP which is the same as the
           value of the smallest subnormal number. Previously we used
           2^-(MANT_DIG - 1) which is too large a value to be useful. Note that we
           can't use ilogb(0), since that isn't a valid thing to do. As a point
           of comparison Java's ulp returns the next normal value e.g.
           2^(1 - MAX_EXP) for ulp(0), but that is not what we want for
           glibc.  */
        /* Fall through...  */
      case FP_SUBNORMAL:
        /* The next closest subnormal value is a constant distance away.  */
        ulp = std::ldexp (1.0, std::numeric_limits<F>::min_exponent
			       - std::numeric_limits<F>::digits);
        break;

      case FP_NORMAL:
        ulp = std::ldexp (1.0, std::ilogb (value)
			  - std::numeric_limits<F>::digits + 1);
        break;

      default:
	std::unreachable();
        break;
    }
  return ulp;
}

/* Returns the number of ulps that GIVEN is away from EXPECTED.  */
#define ULPDIFF(given, expected) \
        (FUNC(fabs) ((given) - (expected)) / ulp (expected))
template <typename F> F ulpdiff (F given, F expected)
{
  return std::fabs ((given - expected) / ulp (expected));
}

typedef std::map<double, uint64_t> ulpacc_t;

static void ulpacc_reduction (ulpacc_t &inout, ulpacc_t &in)
{
  for (const auto &ulp : in)
    inout[ulp.first] += ulp.second;
}

static int
get_max_thread (void)
{
#ifdef _OPENMP
  return omp_get_max_threads ();
#else
  return 1;
#endif
}

static int
get_thread_num (void)
{
#ifdef _OPENMP
  return omp_get_thread_num ();
#else
  return 0;
#endif
}

static inline double parse_range (const std::string &str)
{
  if (str == "pi")
    return std::numbers::pi_v<double>;
  return std::stod (str);
}

int main (int argc, char *argv[])
{
  namespace po = boost::program_options;

  po::options_description desc{"options"};
  desc.add_options()
    ("help,h", "help")
    ("desc,d", po::value<std::string>(), "input description file");

  po::variables_map vm;
  po::store (po::parse_command_line (argc, argv, desc), vm);
  po::notify (vm);

  if (vm.count("help"))
    {
      std::cout << desc << "\n";
      return 1;
    }

  if (!vm.count("desc"))
    {
      std::cerr << "error: specify description file\n";
      return 1;
    }

  std::ifstream file (vm["desc"].as<std::string>());

  boost::property_tree::ptree jsontree;
  boost::property_tree::json_parser::read_json(file, jsontree);

  std::string function = jsontree.get<std::string>("function");
  auto func = get_univariate (function).value();
  auto func_ref = get_univariate_ref (function, FE_TONEAREST).value();

  struct range_t {
    double start;
    double end;
    uint64_t count;
  };
  std::vector<range_t> ranges;
  {
    const boost::property_tree::ptree& ptranges = jsontree.get_child ("ranges");
    for (const auto &ptrange : ptranges)
      {
	double start = parse_range (ptrange.second.get<std::string>("start"));
	double end = parse_range (ptrange.second.get<std::string>("end"));
	uint64_t count = ptrange.second.get<uint64_t>("count");
	ranges.push_back (range_t { start, end, count});
      }
  }

  {
    std::vector<std::mt19937_64> gens (get_max_thread ());
    for (auto &gen : gens)
      gen = seed_rng<std::mt19937_64>();

    ulpacc_t ulpacc;

    #pragma omp declare reduction(ulpacc_reduction :                \
                                  ulpacc_t :                        \
                                  ulpacc_reduction(omp_out, omp_in))\
                    initializer (omp_priv=omp_orig)

    for (const auto &range : ranges)
      {
        std::uniform_real_distribution<double> dist(range.start, range.end);

	ulpacc_t ulpaccrange;

	#pragma omp parallel for shared(dist) reduction(ulpacc_reduction : ulpaccrange)
	for (auto i = 0 ; i < range.count; i++)
	  {
	    float input;

	    input = dist(gens[get_thread_num()]);

	    double computed = func (input);
	    double expected = func_ref (input);

	    double diff = std::fabs (computed - expected);
	    double ulps = ulpdiff (computed, expected);

	    ulpaccrange[ulps] += 1;
	  }

	ulpacc_reduction (ulpacc, ulpaccrange);
      }

    const std::uint64_t ulptotal =
      std::accumulate (ulpacc.begin(),
		       ulpacc.end(),
		       0,
		       [](const uint64_t previous, const std::pair<double, uint64_t> &p)
		       { return previous + p.second; });

    for (const auto &ulp : ulpacc)
      println ("{:g}: {:16} {:6.2f}%",
	       ulp.first,
	       ulp.second,
	       ((double)ulp.second / (double)ulptotal) * 100.0);
  }
}
