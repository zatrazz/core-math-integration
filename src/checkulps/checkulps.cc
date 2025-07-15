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
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include <climits>
#include <fenv.h>
#include <omp.h>

#include "refimpls.h"
#include "wyhash64.h"

using namespace refimpls;
typedef wyhash64 rng_t;

template <typename... Args>
inline void
println (const std::format_string<Args...> fmt, Args &&...args)
{
  std::cout << std::vformat (fmt.get (), std::make_format_args (args...))
            << '\n';
}

template <typename... Args>
inline void
error (const std::format_string<Args...> fmt, Args &&...args)
{
  std::cerr << "error: "
            << std::vformat (fmt.get (), std::make_format_args (args...))
            << '\n';
  std::exit (EXIT_FAILURE);
}

struct round_modes_t
{
  const std::string_view name;
  const std::string_view abbrev;
  int mode;
};

static const std::array k_round_modes = {
#define DEF_ROUND_MODE(__mode, __abbrev)                                      \
  round_modes_t { #__mode, __abbrev, __mode }
  DEF_ROUND_MODE (FE_TONEAREST, "rndn"),
  DEF_ROUND_MODE (FE_UPWARD, "rndu"),
  DEF_ROUND_MODE (FE_DOWNWARD, "rndd"),
  DEF_ROUND_MODE (FE_TOWARDZERO, "rndz"),
#undef DEF_ROUND_MODE
};

static const std::string k_all_round_modes_option = "rndn,rndu,rndd,rndz";

static std::vector<std::string>
split_with_ranges (const std::string_view &s, std::string_view delimiter)
{
  std::vector<std::string> tokens;
  for (const auto &subrange : s | std::views::split (delimiter))
    tokens.push_back (std::string (subrange.begin (), subrange.end ()));
  return tokens;
}

typedef std::vector<round_modes_t> round_set;

static round_set
round_from_option (const std::string_view &rnds)
{
  auto round_modes = split_with_ranges (rnds, ",");
  round_set ret;
  std::copy_if (
      k_round_modes.begin (), k_round_modes.end (), std::back_inserter (ret),
      [&round_modes] (const round_modes_t &r) {
        return std::find_if (
                   round_modes.begin (), round_modes.end (),
                   [&r] (const auto &rnd) { return r.abbrev == rnd; })
               != round_modes.end ();
      });
  return ret;
}

class scope_rouding_t
{
public:
  explicit scope_rouding_t (int rnd)
  {
    curr_rnd = fegetround ();
    assert (fesetround (rnd) == 0);
  }
  ~scope_rouding_t () { assert (fesetround (curr_rnd) == 0); }

  // Prevent copying and moving to ensure single execution
  scope_rouding_t (const scope_rouding_t &) = delete;
  scope_rouding_t &operator= (const scope_rouding_t &) = delete;
  scope_rouding_t (scope_rouding_t &&) = delete;
  scope_rouding_t &operator= (scope_rouding_t &&) = delete;

private:
  int curr_rnd;
};
enum class fail_mode_t
{
  none,
  first,
};

static std::expected<fail_mode_t, bool>
fail_mode_from_options (const std::string_view &modename)
{
  struct fail_mode_desc_t
  {
    std::string name;
    fail_mode_t mode;
  };
  static const std::vector<fail_mode_desc_t> k_fail_modes
      = { { "none", fail_mode_t::none }, { "first", fail_mode_t::first } };
  auto it = std::find_if (k_fail_modes.begin (), k_fail_modes.end (),
                          [&modename] (const fail_mode_desc_t &elem) {
                            return elem.name == modename;
                          });
  if (it != k_fail_modes.end ())
    return (*it).mode;
  return std::unexpected (false);
}

/* Returns the size of an ulp for VALUE.  */
template <typename F>
F
ulp (F value)
{
  F ulp;

  switch (std::fpclassify (value))
    {
    case FP_ZERO:
      /* We compute the distance to the next FP which is the same as the
         value of the smallest subnormal number. Previously we used
         2^-(MANT_DIG - 1) which is too large a value to be useful. Note that
         we can't use ilogb(0), since that isn't a valid thing to do. As a
         point of comparison Java's ulp returns the next normal value e.g. 2^(1
         - MAX_EXP) for ulp(0), but that is not what we want for glibc.  */
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
      std::unreachable ();
      break;
    }
  return ulp;
}

/* Returns the number of ulps that GIVEN is away from EXPECTED.  */
#define ULPDIFF(given, expected)                                              \
  (FUNC (fabs) ((given) - (expected)) / ulp (expected))
template <typename F>
F
ulpdiff (F given, F expected)
{
  return std::fabs (given - expected) / ulp (expected);
}

typedef std::map<double, uint64_t> ulpacc_t;

static void
ulpacc_reduction (ulpacc_t &inout, ulpacc_t &in)
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

static inline double
parse_range (const std::string &str)
{
  if (str == "-pi")
    return -std::numbers::pi_v<double>;
  else if (str == "pi")
    return std::numbers::pi_v<double>;
  else if (str == "2pi")
    return 2.0 * std::numbers::pi_v<double>;
  else if (str == "min")
    return std::numeric_limits<double>::min ();
  else if (str == "-min")
    return -std::numeric_limits<double>::min ();
  else if (str == "max")
    return std::numeric_limits<double>::max ();
  else if (str == "-max")
    return -std::numeric_limits<double>::max ();
  return std::stod (str);
}

struct range_t
{
  double start;
  double end;
  uint64_t count;
};

static void
print_acc (const std::string_view &rndname, const range_t &range,
           const ulpacc_t &ulpacc)
{
  const std::uint64_t ulptotal = std::accumulate (
      ulpacc.begin (), ulpacc.end (), 0,
      [] (const uint64_t previous, const std::pair<double, uint64_t> &p) {
        return previous + p.second;
      });

  println ("Checking rounding mode {:13}, range [{:9.2g},{:9.2g}], count {}",
           rndname, range.start, range.end, ulptotal);

  for (const auto &ulp : ulpacc)
    println ("    {:g}: {:16} {:6.2f}%", ulp.first, ulp.second,
             ((double)ulp.second / (double)ulptotal) * 100.0);
}

static std::vector<rng_t::state_type> rng_states;

static void
init_random_state (void)
{
  std::random_device rd{ "/dev/random" };
  rng_states.resize (get_max_thread ());
  for (auto &s : rng_states)
    // std::random_device max is UINT32_MAX
    s = (rng_t::state_type)rd () << 32 | rd ();
}

template <typename F> struct result_t
{
  typedef F float_type;

  result_t (int r, F c, F e) : rnd (r), computed (c), expected (e)
  {
    ulp = ulpdiff (computed, expected);
  }

  virtual bool
  check (fail_mode_t failmode) const
  {
    if (failmode == fail_mode_t::first && ulp >= 1.0)
      return false;
    return true;
  }

  const std::string_view
  rounding_string () const
  {
    switch (rnd)
      {
      case FE_TONEAREST:
        return "FE_TONEAREST";
      case FE_UPWARD:
        return "FE_UPWARD";
      case FE_DOWNWARD:
        return "FE_DOWNWARD";
      case FE_TOWARDZERO:
        return "FE_TOWARDZERO";
      default:
        std::unreachable ();
      }
  }

  virtual std::ostream &print (std::ostream &) const = 0;

  friend std::ostream &
  operator<< (std::ostream &os, const result_t &e)
  {
    return e.print (os);
  }

  int rnd;
  F computed;
  F expected;
  F ulp;
};

template <typename F> struct result_univariate_t : public result_t<F>
{
public:
  explicit result_univariate_t (int r, F i, F c, F e)
      : result_t<F> (r, c, e), input (i)
  {
  }

  virtual std::ostream &
  print (std::ostream &os) const
  {
    os << std::format (
        "{} ulp={:1.0f} input={:a} computed={:a} expected={:a}\n",
        result_t<F>::rounding_string (), result_t<F>::ulp, input,
        result_t<F>::computed, result_t<F>::expected);
    return os;
  }

  F input;
};

typedef result_univariate_t<float> result_univariate_binary32_t;
typedef result_univariate_t<double> result_univariate_binary64_t;

template <typename F> struct result_bivariate_t : public result_t<F>
{
public:
  explicit result_bivariate_t (int r, F i0, F i1, F c, F e)
      : result_t<F> (r, c, e), input0 (i0), input1 (i1)
  {
  }

  virtual std::ostream &
  print (std::ostream &os) const
  {
    os << std::format (
        "{} ulp={:1.0f} input=({:a},{:a}) computed={:a} expected={:a}\n",
        result_t<F>::rounding_string (), result_t<F>::ulp, input0, input1,
        result_t<F>::computed, result_t<F>::expected);
    return os;
  }

  F input0;
  F input1;
};

typedef result_bivariate_t<float> result_bivariate_binary32_t;
typedef result_bivariate_t<double> result_bivariate_binary64_t;

template <typename RET> struct sample_t
{
  virtual std::unique_ptr<RET>
  operator() (rng_t &, std::uniform_real_distribution<typename RET::float_type> &, int) const
      = 0;
};

template <typename FUNC, typename FUNC_REF, typename RET>
class sample_univariate_t : public sample_t <RET>
{
  typedef typename RET::float_type float_type;

public:
  sample_univariate_t (FUNC &f, FUNC_REF &ref_f) : func (f), ref_func (ref_f)
  {
  }

  std::unique_ptr<RET>
  operator() (rng_t &gen, std::uniform_real_distribution<float_type> &dist,
              int rnd) const
  {
    float_type input = dist (gen);
    float_type computed = func (input);
    float_type expected = ref_func (input, rnd);

    return std::make_unique<RET> (rnd, input, computed, expected);
  }

private:
  const FUNC &func;
  const FUNC_REF &ref_func;
};

typedef sample_univariate_t<univariate_binary32_t, univariate_binary32_ref_t,
                            result_univariate_binary32_t>
    sample_univariate_binary32_t;
typedef sample_univariate_t<univariate_binary64_t, univariate_binary64_ref_t,
                            result_univariate_binary64_t>
    sample_univariate_binary64_t;


template <typename FUNC, typename FUNC_REF, typename RET>
class sample_bivariate_t : public sample_t<RET>
{
  typedef typename RET::float_type float_type;

public:
  sample_bivariate_t (FUNC &f, FUNC_REF &ref_f) : func (f), ref_func (ref_f)
  {
  }

  std::unique_ptr<RET>
  operator() (rng_t &gen, std::uniform_real_distribution<float_type> &dist,
              int rnd) const
  {
    float_type input0 = dist (gen);
    float_type input1 = dist (gen);
    float_type computed = func (input0, input1);
    float_type expected = ref_func (input0, input1, rnd);

    return std::make_unique<RET> (rnd, input0, input1, computed, expected);
  }

private:
  const FUNC &func;
  const FUNC_REF &ref_func;
};

typedef sample_bivariate_t<bivariate_binary32_t, bivariate_binary32_ref_t,
                            result_bivariate_binary32_t>
    sample_bivariate_binary32_t;
typedef sample_bivariate_t<bivariate_binary64_t, bivariate_binary64_ref_t,
                            result_bivariate_binary64_t>
    sample_bivariate_binary64_t;

template <typename RET>
static void
check_variate (const sample_t<RET> &sample, const std::vector<range_t> &ranges,
               const round_set &round_modes, fail_mode_t failmode)
{
  std::vector<rng_t> gens (rng_states.size ());

  for (auto &rnd : round_modes)
    {
      /* Reseed the RNG generator with the same seed to check the same numbers
         for each rounding mode.  */
      for (unsigned i = 0; i < rng_states.size (); i++)
        gens[i] = rng_t (rng_states[i]);

#pragma omp declare reduction(ulpacc_reduction                                \
:ulpacc_t : ulpacc_reduction(omp_out, omp_in))                                \
    initializer(omp_priv = omp_orig)

      for (const auto &range : ranges)
        {
          std::uniform_real_distribution<typename RET::float_type> dist (range.start, range.end);

          ulpacc_t ulpaccrange;

#pragma omp parallel firstprivate(dist, failmode) shared(sample, rnd)
          {
            scope_rouding_t scope_rounding (rnd.mode);

#pragma omp for reduction(ulpacc_reduction : ulpaccrange)
            for (unsigned i = 0; i < range.count; i++)
              {
                auto ret = sample (gens[get_thread_num ()], dist, rnd.mode);
                if (!ret->check (failmode))
#pragma omp critical
                  {
                    std::cerr << *ret;
                    std::exit (EXIT_FAILURE);
                  }
                ulpaccrange[ret->ulp] += 1;
              }
          }

          print_acc (rnd.name, range, ulpaccrange);
        }

      println ("");
    }
}

int
main (int argc, char *argv[])
{
  namespace po = boost::program_options;

  po::options_description desc{ "options" };
  desc.add_options () ("help,h", "help") (
      "verbose,v", po::bool_switch ()->default_value (false)) (
      "core,c", po::bool_switch ()->default_value (false),
      "check CORE-MATH routines") ("desc,d", po::value<std::string> (),
                                   "input description file") (
      "rnd,r",
      po::value<std::string> ()->default_value (k_all_round_modes_option),
      "rounding mode to test") (
      "fail,f", po::value<std::string> ()->default_value ("none"));

  po::variables_map vm;
  po::store (po::parse_command_line (argc, argv, desc), vm);
  po::notify (vm);

  if (vm.count ("help"))
    {
      std::cout << desc << "\n";
      return 1;
    }

  if (!vm.count ("desc"))
    error ("no description file\n");

  bool verbose = vm["verbose"].as<bool> ();
  bool coremath = vm["core"].as<bool> ();

  const round_set round_modes
      = round_from_option (vm["rnd"].as<std::string> ());

  fail_mode_t failmode
      = fail_mode_from_options (vm["fail"].as<std::string> ()).value ();

  std::ifstream file (vm["desc"].as<std::string> ());

  boost::property_tree::ptree jsontree;
  boost::property_tree::json_parser::read_json (file, jsontree);

  std::string function = jsontree.get<std::string> ("function");

  std::vector<range_t> ranges;
  {
    const boost::property_tree::ptree &ptranges
        = jsontree.get_child ("ranges");
    for (const auto &ptrange : ptranges)
      {
        double start = parse_range (ptrange.second.get<std::string> ("start"));
        double end = parse_range (ptrange.second.get<std::string> ("end"));
        uint64_t count = ptrange.second.get<uint64_t> ("count");
        ranges.push_back (range_t{ start, end, count });
        if (verbose)
          println ("range=[start={:a},end={:a},count={}", start, end, count);
      }
  }

  auto type = get_func_type (function);
  if (!type)
    error ("invalid function: {}", function);

  init_random_state ();

  switch (type.value ())
    {
    case refimpls::func_type_t::binary32_univariate:
      {
        auto func = get_univariate_binary32 (function, coremath).value ();
        if (!func.first)
          error ("libc does not provide {}", function);

        sample_univariate_binary32_t sample{ func.first, func.second };
        check_variate (sample, ranges, round_modes, failmode);
      }
      break;
    case refimpls::func_type_t::binary32_bivariate:
      {
        auto func = get_bivariate_binary32 (function, coremath).value ();
        if (!func.first)
          error ("libc does not provide {}", function);

        sample_bivariate_binary32_t sample{ func.first, func.second };
        check_variate (sample, ranges, round_modes, failmode);
      }
      break;
    case refimpls::func_type_t::binary64_univariate:
      {
        auto func = get_univariate_binary64 (function, coremath).value ();
        if (!func.first)
          error ("libc does not provide {}", function);

        sample_univariate_binary64_t sample{ func.first, func.second };
        check_variate (sample, ranges, round_modes, failmode);
      }
      break;
    case refimpls::func_type_t::binary64_bivariate:
      {
        auto func = get_bivariate_binary64 (function, coremath).value ();
        if (!func.first)
          error ("libc does not provide {}", function);

        sample_bivariate_binary64_t sample{ func.first, func.second };
        check_variate (sample, ranges, round_modes, failmode);
      }
      break;
    default:
      error ("function {} not implemented", function);
    }
}
