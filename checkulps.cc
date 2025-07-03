#include <algorithm>
#include <cassert>
#include <format>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <numbers>
#include <random>
#include <ranges>
#include <type_traits>

#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include <climits>
#include <fenv.h>
#include <omp.h>

#include "wyhash64.h"
#include "refimpls.h"

using namespace refimpls;
typedef wyhash64 rng_t;

template<typename... Args>
inline void println(const std::format_string<Args...> fmt, Args&&... args)
{
  std::cout << std::vformat(fmt.get(), std::make_format_args(args...)) << '\n';
}

template<typename... Args>
inline void error(const std::format_string<Args...> fmt, Args&&... args)
{
  std::cerr << "error: " << std::vformat(fmt.get(), std::make_format_args(args...)) << '\n';
  std::exit (EXIT_FAILURE);
}

struct round_modes_t
{
  std::string name;
  std::string abbrev;
  int mode;
};
typedef std::vector<round_modes_t> round_set;

static const round_set k_round_modes =
{
#define DEF_ROUND_MODE(__mode, __abbrev) { #__mode, __abbrev, __mode }
  DEF_ROUND_MODE (FE_TONEAREST,  "rndn"),
  DEF_ROUND_MODE (FE_UPWARD,     "rndu"),
  DEF_ROUND_MODE (FE_DOWNWARD,   "rndd"),
  DEF_ROUND_MODE (FE_TOWARDZERO, "rndz"),
#undef DEF_ROUND_MODE
};


static std::vector<std::string_view>
split_with_ranges(const std::string_view& s, std::string_view delimiter)
{
  std::vector<std::string_view> tokens;
  for (const auto& subrange : s | std::views::split(delimiter))
    tokens.push_back(std::string(subrange.begin(), subrange.end()));
  return tokens;
}

static round_set round_from_option (const std::string_view& rnds)
{
  if (rnds == "all")
    return k_round_modes;

  auto round_modes = split_with_ranges (rnds, ",");
  round_set ret;
  std::copy_if(k_round_modes.begin(),
	       k_round_modes.end(),
	       std::back_inserter(ret),
	       [&round_modes](const round_modes_t &r) {
	         return std::find (round_modes.begin(),
				   round_modes.end(),
				   r.abbrev) != round_modes.end(); });
  return ret;
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

class ScopeRounding {
public:
  explicit ScopeRounding(int rnd) {
    curr_rnd = fegetround ();
    assert (fesetround (rnd) == 0);
  }
  ~ScopeRounding() {
    assert (fesetround (curr_rnd) == 0);
  }

  // Prevent copying and moving to ensure single execution
  ScopeRounding (const ScopeRounding &) = delete;
  ScopeRounding &operator=(const ScopeRounding &) = delete;
  ScopeRounding (ScopeRounding &&) = delete;
  ScopeRounding &operator=(ScopeRounding &&) = delete;

private:
  int curr_rnd;
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


  println ("Checking rounding mode {:13}, range [{:9.2g},{:9.2g}], count {}",
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

static std::vector<rng_t::state_type> rng_states;

static void
init_random_state (void)
{
  std::random_device rd {"/dev/random"};
  rng_states.resize(get_max_thread());
  for (auto &s : rng_states)
    s = (rng_t::state_type)rd() << 32 | rd ();
}

static void
check_univariate (univariate_t func,
		  univariate_ref_t func_ref,
		  const std::vector<range_t> &ranges,
		  const round_set& round_modes)
{
  std::vector<rng_t> gens(rng_states.size());

  for (auto &rnd : round_modes)
    {
      /* Reseed the RNG generator with the same seed to check the same numbers
	 for each rounding mode.  */
      for (unsigned i = 0; i < rng_states.size(); i++)
	gens[i] = rng_t(rng_states[i]);

      ScopeRounding scope_rnd (rnd.mode);

      #pragma omp declare reduction(ulpacc_reduction :                \
				    ulpacc_t :                        \
				    ulpacc_reduction(omp_out, omp_in))\
		      initializer (omp_priv=omp_orig)

      for (const auto &range : ranges)
	{
	  std::uniform_real_distribution<double> dist (range.start, range.end);

	  ulpacc_t ulpaccrange;

	  #pragma omp parallel for reduction(ulpacc_reduction : ulpaccrange) \
				   firstprivate (rnd, dist)
	  for (unsigned i = 0 ; i < range.count; i++)
	    {
	      if (i == 0)
		println ("rnd=({},{},{})\n", rnd.name, rnd.abbrev, rnd.mode);
	      double input = dist(gens[get_thread_num()]);

	      double computed = func (input);
	      double expected = func_ref (input, rnd.mode);

	      double ulps = ulpdiff (computed, expected);

	      ulpaccrange[ulps] += 1;
	    }

	  print_acc (rnd.name, range, ulpaccrange);
	}

      println ("");
    }
}

static void
check_bivariate (bivariate_t func,
		 bivariate_ref_t func_ref,
		 const std::vector<range_t> &ranges,
		 const round_set& round_modes)
{
  std::vector<rng_t> gens(rng_states.size());

  for (auto &rnd : round_modes)
    {
      for (unsigned i = 0; i < rng_states.size(); i++)
	gens[i] = rng_t(rng_states[i]);

      ScopeRounding scope_rnd (rnd.mode);

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
	  for (unsigned i = 0 ; i < range.count; i++)
	    {
	      double input0 = dist(gens[get_thread_num()]);
	      double input1 = dist(gens[get_thread_num()]);

	      double computed = func (input0, input1);
	      double expected = func_ref (input0, input1, rnd.mode);

	      double ulps = ulpdiff (computed, expected);

	      ulpaccrange[ulps] += 1;
	    }

	  print_acc (rnd.name, range, ulpaccrange);
	}

      println ("");
    }
}

int main (int argc, char *argv[])
{
  namespace po = boost::program_options;

  po::options_description desc{"options"};
  desc.add_options()
    ("help,h",    "help")
    ("verbose,v", po::bool_switch()->default_value(false))
    ("core,c",    po::bool_switch()->default_value(false), "check CORE-MATH routines")
    ("desc,d",    po::value<std::string>(), "input description file")
    ("rnd,r",     po::value<std::string>()->default_value("all"), "rounding mode to test");

  po::variables_map vm;
  po::store (po::parse_command_line (argc, argv, desc), vm);
  po::notify (vm);

  if (vm.count("help"))
    {
      std::cout << desc << "\n";
      return 1;
    }

  if (!vm.count("desc"))
    error ("no description file\n");

  bool verbose = vm["verbose"].as<bool>();
  bool coremath = vm["core"].as<bool>();

  const round_set round_modes = round_from_option (vm["rnd"].as<std::string>());

  std::ifstream file (vm["desc"].as<std::string>());

  boost::property_tree::ptree jsontree;
  boost::property_tree::json_parser::read_json(file, jsontree);

  std::string function = jsontree.get<std::string>("function");

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

  auto type = get_func_type (function);
  if (!type)
    error ("invalid function: {}", function);

  init_random_state ();

  switch (type.value())
  {
  case refimpls::func_type_t::univariate:
    {
      auto func = get_univariate (function, coremath).value();
      if (!func.first)
	error ("function {} not provided by libc\n", function);

      check_univariate (func.first, func.second, ranges, round_modes);
    } break;
  case refimpls::func_type_t::bivariate:
    {
      auto func = get_bivariate (function, coremath).value();
      if (!func.first)
	error ("function {} not provided by libc\n", function);

      check_bivariate (func.first, func.second, ranges, round_modes);
    } break;
  }
}
