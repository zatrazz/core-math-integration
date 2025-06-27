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
  else if (str == "min")
    return std::numeric_limits<double>::min();
  else if (str == "-min")
    return -std::numeric_limits<double>::min();
  else if (str == "max")
    return std::numeric_limits<double>::max();
  else if (str == "-max")
    return -std::numeric_limits<double>::max();
  return std::stod (str);
}

struct range_t {
  double start;
  double end;
  uint64_t count;
};

struct round_modes_t {
  std::string name;
  int mode;
};

static void
print_acc (const std::string &rndname,
	   const range_t& range,
	   const ulpacc_t& ulpacc)
{
  const std::uint64_t ulptotal =
    std::accumulate (ulpacc.begin(),
		     ulpacc.end(),
		     0,
		     [](const uint64_t previous, const std::pair<double, uint64_t> &p)
		     { return previous + p.second; });


  println ("Checking rounding mode {:13}, range [{:8.2g},{:8.2g}], count {}",
	   rndname,
	   range.start,
	   range.end,
	   range.count);

  for (const auto &ulp : ulpacc)
    println ("    {:g}: {:16} {:6.2f}%",
	     ulp.first,
	     ulp.second,
	     ((double)ulp.second / (double)ulptotal) * 100.0);
}

static void
check_univariate (univariate_t func,
		  univariate_t func_ref,
		  const std::vector<range_t> &ranges,
		  round_modes_t rnd)
{
  std::vector<std::mt19937_64> gens (get_max_thread ());
  for (auto &gen : gens)
    gen = seed_rng<std::mt19937_64>();

  if (fesetround (rnd.mode) != 0)
    {
      std::cerr << "error: invalid rounding mode: " << rnd.mode;
      std::abort ();
    }

  #pragma omp declare reduction(ulpacc_reduction :                \
				ulpacc_t :                        \
				ulpacc_reduction(omp_out, omp_in))\
		  initializer (omp_priv=omp_orig)

  for (const auto &range : ranges)
    {
      ulpacc_t ulpaccrange;

      #pragma omp paralle
      {
        std::uniform_real_distribution<double> dist (range.start, range.end);

	#pragma omp parallel for reduction(ulpacc_reduction : ulpaccrange)
	for (auto i = 0 ; i < range.count; i++)
	  {
	    double input = input = dist(gens[get_thread_num()]);

	    double computed = func (input);
	    double expected = func_ref (input);

#if 0
	    if (!std::isnormal (computed) || !std::isnormal (expected))
	      {
		println ("input={0:a},{0:g} computed={1:a},{1:g} "
			 "expected={2:a},{2:g}\n",
			 input,
			 computed,
			 expected);
		std::exit (1);
	      }
#endif

	    double diff = std::fabs (computed - expected);
	    double ulps = ulpdiff (computed, expected);

	    ulpaccrange[ulps] += 1;
	  }

      }

      print_acc (rnd.name, range, ulpaccrange);
    }

  println ("");
}

int main (int argc, char *argv[])
{
  namespace po = boost::program_options;

  po::options_description desc{"options"};
  desc.add_options()
    ("help,h",    "help")
    ("verbose,v", po::bool_switch()->default_value(false))
    ("desc,d",    po::value<std::string>(), "input description file");

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

  bool verbose = vm["verbose"].as<bool>();

  std::ifstream file (vm["desc"].as<std::string>());

  boost::property_tree::ptree jsontree;
  boost::property_tree::json_parser::read_json(file, jsontree);

  std::string function = jsontree.get<std::string>("function");
  auto func = get_univariate (function).value();
  auto func_ref = get_univariate_ref (function, FE_TONEAREST).value();

  std::vector<range_t> ranges;
  {
    const boost::property_tree::ptree& ptranges = jsontree.get_child ("ranges");
    for (const auto &ptrange : ptranges)
      {
	double start = parse_range (ptrange.second.get<std::string>("start"));
	double end = parse_range (ptrange.second.get<std::string>("end"));
	uint64_t count = ptrange.second.get<uint64_t>("count");
	ranges.push_back (range_t { start, end, count});
	if (verbose)
	  println ("range=[start={:a},end={:a},count={}",
		   start, end, count);
      }
  }

  std::vector<round_modes_t> round_modes = {
    { "FE_TONEAREST",  FE_TONEAREST },
    { "FE_UPWARD"   ,  FE_TONEAREST },
    { "FE_DOWNWARD",   FE_DOWNWARD },
    { "FE_TOWARDZERO", FE_TOWARDZERO },
  };

  for (auto &rnd : round_modes)
    check_univariate (func, func_ref, ranges, rnd);
}
