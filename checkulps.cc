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

#include "pcg_random.hpp"
#include "refimpls.h"

using namespace refimpls;
typedef pcg64 rng_t;
typedef pcg_extras::seed_seq_from<std::random_device> rng_seed_t;

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

class ScopeGuard {
public:
  explicit ScopeGuard(std::function<void()> f)
      : onexit(std::move(f)) {}
  ~ScopeGuard() {
    onexit();
  }

  // Prevent copying and moving to ensure single execution
  ScopeGuard (const ScopeGuard &) = delete;
  ScopeGuard &operator=(const ScopeGuard &) = delete;
  ScopeGuard (ScopeGuard &&) = delete;
  ScopeGuard &operator=(ScopeGuard &&) = delete;

private:
  std::function<void()> onexit;
};

static inline double parse_range (const std::string &str)
{
  if (str == "-pi")
    return -std::numbers::pi_v<double>;
  else if (str == "pi")
    return std::numbers::pi_v<double>;
  else if (str == "2pi")
    return 2.0 * std::numbers::pi_v<double>;
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
		  univariate_ref_t func_ref,
		  const std::vector<range_t> &ranges,
		  const round_modes_t rnd)
{
  int current_rnd = fesetround (rnd.mode);
  if (current_rnd != 0)
    {
      std::cerr << "error: invalid rounding mode: " << rnd.mode;
      std::abort ();
    }
  ScopeGuard cleanup_rnd ([&current_rnd](){ fesetround (current_rnd); });

  std::vector<rng_t> gens(get_max_thread());
  for (auto &g : gens)
    g = rng_t(rng_seed_t{});

  #pragma omp declare reduction(ulpacc_reduction :                \
				ulpacc_t :                        \
				ulpacc_reduction(omp_out, omp_in))\
		  initializer (omp_priv=omp_orig)

  for (const auto &range : ranges)
    {
      std::uniform_real_distribution<double> dist (range.start, range.end);

      ulpacc_t ulpaccrange;

      #pragma omp parallel for reduction(ulpacc_reduction : ulpaccrange) \
			       firstprivate (dist)
      for (auto i = 0 ; i < range.count; i++)
	{
	  double input = input = dist(gens[get_thread_num()]);

	  double computed = func (input);
	  double expected = func_ref (input, rnd.mode);

	  double diff = std::fabs (computed - expected);
	  double ulps = ulpdiff (computed, expected);

	  ulpaccrange[ulps] += 1;
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
    ("core,c",    po::bool_switch()->default_value(false), "check CORE-MATH routines")
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
  bool coremath = vm["core"].as<bool>();

  std::ifstream file (vm["desc"].as<std::string>());

  boost::property_tree::ptree jsontree;
  boost::property_tree::json_parser::read_json(file, jsontree);

  std::string function = jsontree.get<std::string>("function");

  auto func = get_univariate (function, coremath).value();
  if (!func)
    {
      std::cerr << std::format("error: function {} not provided by glibc\n",
			       function);
      return 1;
    }
  auto func_ref = get_univariate_ref (function).value();

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

  static const std::vector<round_modes_t> round_modes = {
    { "FE_TONEAREST",  FE_TONEAREST },
    { "FE_UPWARD"   ,  FE_TONEAREST },
    { "FE_DOWNWARD",   FE_DOWNWARD },
    { "FE_TOWARDZERO", FE_TOWARDZERO },
  };

  for (auto &rnd : round_modes)
    check_univariate (func, func_ref, ranges, rnd);
}
