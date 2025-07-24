#include <algorithm>
#include <numbers>
#include <random>
#include <ranges>
#include <chrono>
#include <iostream>

#include <boost/program_options.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

#include <fenv.h>
#include <omp.h>

#include "refimpls.h"
#include "wyhash64.h"

using namespace refimpls;
typedef wyhash64 rng_t;

using clock_type = std::chrono::high_resolution_clock;
using duration_type = std::chrono::duration<double, std::milli>;

template<typename... Args>
void println_ts (std::format_string<Args...> fmt, Args&&... args)
{
  auto now = std::chrono::system_clock::now();
  auto seconds = std::chrono::floor<std::chrono::seconds>(now);
  std::print("[{:%Y-%m-%d %H:%M:%S}] ", seconds);
  std::println(fmt, std::forward<Args>(args)...);
}

// Information class used to generate full ranges, mainly for testing
// all binary32 normal and subnormal numbers.
template <typename T> struct float_ranges_t
{
  static constexpr uint64_t plus_normal_min = 0;
  static constexpr uint64_t plus_normal_max = 0;
  static constexpr uint64_t plus_subnormal_min = 0;
  static constexpr uint64_t plus_subnormal_max = 0;
  static constexpr uint64_t neg_normal_min = 0;
  static constexpr uint64_t neg_normal_max = 0;
  static constexpr uint64_t neg_subnormal_min = 0;
  static constexpr uint64_t neg_subnormal_max = 0;
};

template <> struct float_ranges_t<float>
{
  static constexpr uint64_t plus_normal_min = 0x00800000;
  static constexpr uint64_t plus_normal_max = 0x7F7FFFFF;
  static constexpr uint64_t plus_subnormal_min = 0x00000001;
  static constexpr uint64_t plus_subnormal_max = 0x007FFFFF;
  static constexpr uint64_t neg_normal_min = 0x80800000;
  static constexpr uint64_t neg_normal_max = 0xFF7FFFFF;
  static constexpr uint64_t neg_subnormal_min = 0x80000001;
  static constexpr uint64_t neg_subnormal_max = 0x807FFFFF;

  static constexpr float
  from_integer (uint32_t u)
  {
    union
    {
      uint32_t u;
      float f;
    } r = { .u = u };
    return r.f;
  }
};

template <> struct float_ranges_t<double>
{
  static constexpr uint64_t plus_normal_min = UINT64_C(0x0008000000000000);
  static constexpr uint64_t plus_normal_max = UINT64_C(0x7FEFFFFFFFFFFFFF);
  static constexpr uint64_t plus_subnormal_min = UINT64_C(0x0000000000000001);
  static constexpr uint64_t plus_subnormal_max = UINT64_C(0x000FFFFFFFFFFFFF);
  static constexpr uint64_t neg_normal_min = UINT64_C(0x8008000000000000);
  static constexpr uint64_t neg_normal_max = UINT64_C(0xFFEFFFFFFFFFFFFF);
  static constexpr uint64_t neg_subnormal_min = UINT64_C(0x8000000000000001);
  static constexpr uint64_t neg_subnormal_max = UINT64_C(0x800FFFFFFFFFFFFF);

  static constexpr double
  from_integer (uint64_t u)
  {
    union
    {
      uint64_t u;
      double f;
    } r = { .u = u };
    return r.f;
  }
};

template <std::floating_point T> bool issignaling (T);

template <>
bool
issignaling<float> (float x)
{
  union
  {
    float f;
    std::uint32_t u;
  } n = { .f = x };
  n.u ^= 0x00400000;
  return (n.u & 0x7fffffff) > 0x7fc00000;
}

template <>
bool
issignaling<double> (double x)
{
  union
  {
    double f;
    std::uint64_t u;
  } n = { .f = x };
  n.u ^= UINT64_C (0x0008000000000000);
  return (n.u & UINT64_C (0x7fffffffffffffff)) > UINT64_C (0x7ff8000000000000);
}

template <typename... Args>
[[noreturn]] inline void
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
  int mode;

  bool
  operator== (const round_modes_t &other) const
  {
    return mode == other.mode;
  }
};


typedef std::vector<round_modes_t> round_set;

static const std::map<std::string_view, round_modes_t> k_round_modes = {
#define DEF_ROUND_MODE(__mode, __abbrev)                                      \
  {                                                                           \
    __abbrev, round_modes_t { #__mode, __mode }                               \
  }
  DEF_ROUND_MODE (FE_TONEAREST, "rndn"),
  DEF_ROUND_MODE (FE_UPWARD, "rndu"),
  DEF_ROUND_MODE (FE_DOWNWARD, "rndd"),
  DEF_ROUND_MODE (FE_TOWARDZERO, "rndz"),
#undef DEF_ROUND_MODE
};

static const std::string k_all_round_modes_option = std::accumulate (
    std::next (k_round_modes.begin ()), k_round_modes.end (),
    k_round_modes.empty () ? "" : std::string (k_round_modes.begin ()->first),
    [] (const std::string &acc, const auto &pair) {
      return acc + "," + std::string (pair.first);
    });

static const std::string_view
round_name (int rnd)
{
  if (auto it
      = std::find_if (k_round_modes.begin (), k_round_modes.end (),
                      [&rnd] (const auto &r) { return rnd == r.second.mode; });
      it != k_round_modes.end ())
    return it->second.name;
  std::unreachable ();
}

static std::vector<std::string>
split_with_ranges (const std::string_view &s, std::string_view delimiter)
{
  std::vector<std::string> tokens;
  for (const auto &subrange : s | std::views::split (delimiter))
    tokens.push_back (std::string (subrange.begin (), subrange.end ()));
  return tokens;
}

static round_set
round_from_option (const std::string_view &rnds)
{
  auto round_modes = split_with_ranges (rnds, ",");
  round_set ret;
  for (auto &rnd : round_modes)
    if (auto it = k_round_modes.find (rnd); it != k_round_modes.end ())
      {
        if (std::ranges::contains (ret, it->second))
          error ("rounding mode already defined: {}", rnd);
        ret.push_back (it->second);
      }
    else
      error ("invalid rounding mode: {}", rnd);
  return ret;
}

template <typename F> class round_setup_t
{
  int saved_rnd;
  std::pair<mpfr_exp_t, mpfr_exp_t> saved_ref_param;

  void
  setup_rnd (int rnd)
  {
    if (fesetround (rnd) != 0)
      error ("fesetround ({}) failed (errno={})", round_name (rnd), errno);
  }

public:
  explicit round_setup_t (int rnd)
  {
    saved_rnd = fegetround ();
    setup_rnd (rnd);
    refimpls::setup_ref_impl<F> ();
  }
  ~round_setup_t () { setup_rnd (saved_rnd); }

  // Prevent copying and moving to ensure single execution
  round_setup_t (const round_setup_t &) = delete;
  round_setup_t &operator= (const round_setup_t &) = delete;
  round_setup_t (round_setup_t &&) = delete;
  round_setup_t &operator= (round_setup_t &&) = delete;
};

enum class fail_mode_t
{
  none,
  first,
};

static const std::map<std::string_view, fail_mode_t> k_fail_modes
    = { { "none", fail_mode_t::none }, { "first", fail_mode_t::first } };

fail_mode_t
fail_mode_from_options (const std::string_view &failmode)
{
  if (auto it = k_fail_modes.find (failmode); it != k_fail_modes.end ())
    return it->second;

  error ("invalid fail mode: {}", failmode);
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
      /* Fall through...  */
    case FP_SUBNORMAL:
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
template <typename F>
F
ulpdiff (F given, F expected)
{
  return std::fabs (given - expected) / ulp (expected);
}

template <typename F> using ulpacc_t = std::map<F, uint64_t>;

template <typename F>
static void
ulpacc_reduction (ulpacc_t<F> &inout, ulpacc_t<F> &in)
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

template <typename F>
static inline F
parse_range (const std::string &str)
{
  if (str == "-pi")
    return -std::numbers::pi_v<F>;
  else if (str == "pi")
    return std::numbers::pi_v<F>;
  else if (str == "2pi")
    return 2.0 * std::numbers::pi_v<F>;
  else if (str == "min")
    return std::numeric_limits<F>::min ();
  else if (str == "-min")
    return -std::numeric_limits<F>::min ();
  else if (str == "max")
    return std::numeric_limits<F>::max ();
  else if (str == "-max")
    return -std::numeric_limits<F>::max ();
  return std::stod (str);
}

// Ranges used to represent either the range list for random sampling or
// all number for a defined class.

template <typename F> struct range_random_t
{
  F start;
  F end;
  uint64_t count;
};
template <typename F>
using range_random_list_t = std::vector<range_random_t<F> >;

struct range_full_t
{
  std::string name;
  uint64_t start;
  uint64_t end;
};
typedef std::vector<range_full_t> range_full_list_t;

template <typename F>
static void
print_acc (const std::string_view &rndname, const range_random_t<F> &range,
           const ulpacc_t<F> &ulpacc)
{
  const std::uint64_t ulptotal = std::accumulate (
      ulpacc.begin (), ulpacc.end (), 0,
      [] (const uint64_t previous, const std::pair<double, uint64_t> &p) {
        return previous + p.second;
      });

  println_ts (
      "Checking rounding mode {:13}, range [{:9.2g},{:9.2g}], count {}",
      rndname, range.start, range.end, ulptotal);

  for (const auto &ulp : ulpacc)
    println_ts ("    {:g}: {:16} {:6.2f}%", ulp.first, ulp.second,
                  ((double)ulp.second / (double)ulptotal) * 100.0);
}

template <typename F>
static void
print_acc (const std::string_view &rndname, const range_full_t &range,
           const ulpacc_t<F> &ulpacc)
{
  const std::uint64_t ulptotal = std::accumulate (
      ulpacc.begin (), ulpacc.end (), 0,
      [] (const uint64_t previous, const std::pair<double, uint64_t> &p) {
        return previous + p.second;
      });

  println_ts ("Checking rounding mode {:13}, {}", rndname, range.name);

  for (const auto &ulp : ulpacc)
    println_ts ("    {:g}: {:16} {:6.2f}%", ulp.first, ulp.second,
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
    if (std::isnan (ulp) || std::isinf (ulp))
      // Do not signal an error if the expected value is NaN/Inf.
      ulp = 0.0;
  }

  virtual bool
  check (fail_mode_t failmode) const
  {
    if (failmode == fail_mode_t::first && ulp >= 1.0)
      return false;
    return true;
  }

  virtual bool
  check_full (fail_mode_t failmode) const
  {
    if (issignaling (computed) || issignaling (expected))
      return failmode == fail_mode_t::first ? false : true;
    else if (std::isnan (computed) && std::isnan (expected))
      return true;
    else if (std::isinf (computed) && std::isinf (expected))
      /* Test for sign of infinities.  */
      return std::signbit (computed) == std::signbit (expected);
    else if (std::isinf (computed) || std::isnan (computed)
             || std::isinf (expected) || std::isnan (expected))
      return failmode == fail_mode_t::first ? false : true;

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
        "{} ulp={:1.0f} input=0x{:a} computed=0x{:a} expected=0x{:a}\n",
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
    os << std::format ("{} ulp={:1.0f} input=(0x{:a},0x{:a}) computed=0x{:a} "
                       "expected=0x{:a}\n",
                       result_t<F>::rounding_string (), result_t<F>::ulp,
                       input0, input1, result_t<F>::computed,
                       result_t<F>::expected);
    return os;
  }

  F input0;
  F input1;
};

typedef result_bivariate_t<float> result_bivariate_binary32_t;
typedef result_bivariate_t<double> result_bivariate_binary64_t;

template <typename RET> struct sample_random_t
{
  virtual std::unique_ptr<RET>
  operator() (rng_t &,
              std::uniform_real_distribution<typename RET::float_type> &,
              int) const
      = 0;
};

template <typename FUNC, typename FUNC_REF, typename RET>
class sample_random_univariate_t : public sample_random_t<RET>
{
  typedef typename RET::float_type float_type;

  const FUNC &func;
  const FUNC_REF &ref_func;

public:
  sample_random_univariate_t (FUNC &f, FUNC_REF &ref_f)
      : func (f), ref_func (ref_f)
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
};

template <typename F>
using random_univariate_t
    = sample_random_univariate_t<univariate_t<F>, univariate_ref_t<F>,
                                 result_univariate_t<F> >;

template <typename FUNC, typename FUNC_REF, typename RET>
class sample_random_bivariate_t : public sample_random_t<RET>
{
  typedef typename RET::float_type float_type;

  const FUNC &func;
  const FUNC_REF &ref_func;

public:
  sample_random_bivariate_t (FUNC &f, FUNC_REF &ref_f)
      : func (f), ref_func (ref_f)
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
};

template <typename F>
using random_bivariate_t
    = sample_random_bivariate_t<bivariate_t<F>, bivariate_ref_t<F>,
                                result_bivariate_t<F> >;

template <typename RET> struct sample_full_t
{
  virtual std::unique_ptr<RET> operator() (uint64_t, int) const = 0;
};

template <typename FUNC, typename FUNC_REF, typename RET>
class sample_full_univariate_t : public sample_full_t<RET>
{
  typedef typename RET::float_type float_type;

  const FUNC &func;
  const FUNC_REF &ref_func;

public:
  sample_full_univariate_t (FUNC &f, FUNC_REF &ref_f)
      : func (f), ref_func (ref_f)
  {
  }

  std::unique_ptr<RET>
  operator() (uint64_t n, int rnd) const
  {
    float_type input = float_ranges_t<float_type>::from_integer (n);
    float_type computed = func (input);
    float_type expected = ref_func (input, rnd);

    return std::make_unique<RET> (rnd, input, computed, expected);
  }
};

template <typename F>
using full_univariate_t
    = sample_full_univariate_t<univariate_t<F>, univariate_ref_t<F>,
                               result_univariate_t<F> >;

static void
print_start (const std::string_view &funcname)
{
  println ("Checking function {}\n", funcname);
}

template <typename RET>
static void
check_random_variate (
    const std::string_view &funcname, const sample_random_t<RET> &sample,
    const range_random_list_t<typename RET::float_type> &ranges,
    const round_set &round_modes, fail_mode_t failmode)
{
  using float_type = typename RET::float_type;

  std::vector<rng_t> gens (rng_states.size ());

  print_start (funcname);

  duration_type duration_total {0};

  for (auto &rnd : round_modes)
    {
      /* Reseed the RNG generator with the same seed to check the same numbers
         for each rounding mode.  */
      for (unsigned i = 0; i < rng_states.size (); i++)
        gens[i] = rng_t (rng_states[i]);

#pragma omp declare reduction(                                                \
        ulpacc_reduction : ulpacc_t<float_type> : ulpacc_reduction(           \
                omp_out, omp_in)) initializer(omp_priv = omp_orig)

      auto start = clock_type::now();

      for (const auto &range : ranges)
        {
          std::uniform_real_distribution<float_type> dist (range.start,
                                                           range.end);

          ulpacc_t<float_type> ulpaccrange;

#pragma omp parallel firstprivate(dist, failmode) shared(sample, rnd)
          {
            round_setup_t<float_type> round_setup (rnd.mode);

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

      auto end = clock_type::now();

      duration_type duration { end - start };
      duration_total += duration;

      println_ts ("Elapsed time {}",
		  std::chrono::duration_cast<std::chrono::duration<double>> (duration));

      println_ts ("");
    }


  println_ts ("Total elapsed time {}",
	      std::chrono::duration_cast<std::chrono::duration<double>>(duration_total));
}

template <typename RET>
static void
check_full_variate (const std::string_view &funcname,
                    const sample_full_t<RET> &sample,
                    const range_full_list_t &ranges,
                    const round_set &round_modes, fail_mode_t failmode)
{
  using float_type = typename RET::float_type;

  print_start (funcname);

  for (auto &rnd : round_modes)
    {
#pragma omp declare reduction(                                                \
        ulpacc_reduction : ulpacc_t<float_type> : ulpacc_reduction(           \
                omp_out, omp_in)) initializer(omp_priv = omp_orig)

      for (const auto &range : ranges)
        {
          ulpacc_t<float_type> ulpaccrange;

#pragma omp parallel firstprivate(failmode) shared(sample, rnd)
          {
            round_setup_t<float_type> round_setup (rnd.mode);

// Out of range inputs might take way less time than normal one, also use
// a large chunk size to minimize the overhead from dynamic scheduline.
#pragma omp for reduction(ulpacc_reduction : ulpaccrange) schedule(dynamic)
            for (uint64_t i = range.start; i < range.end; i++)
              {
                auto ret = sample (i, rnd.mode);
                if (!ret->check_full (failmode))
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

      println_ts ("");
    }
}

template <typename F>
using ranges_types_t = std::variant<range_random_list_t<F>, range_full_list_t>;
enum range_type_mode_t
{
  RANGE_RANDOM = 0,
  RANGE_FULL = 1
};

template <typename F>
static ranges_types_t<F>
parse_ranges (const boost::property_tree::ptree &jsontree, bool verbose)
{
  range_random_list_t<F> ranges;
  const boost::property_tree::ptree &ptranges = jsontree.get_child ("ranges");
  for (const auto &ptrange : ptranges)
    {
      F start = parse_range<F> (ptrange.second.get<std::string> ("start"));
      F end = parse_range<F> (ptrange.second.get<std::string> ("end"));
      uint64_t count = ptrange.second.get<uint64_t> ("count");
      ranges.push_back (range_random_t{ start, end, count });
      if (verbose)
        println_ts ("range=[start={:a},end={:a},count={}", start, end,
                      count);
    }
  if (ptranges.size () != 0)
    return ranges;

  range_full_list_t fullranges;
  auto ptfullranges = split_with_ranges (ptranges.get<std::string> (""), ",");
  for (const auto &ptrange : ptfullranges)
    {
      if (ptrange == "normal")
        {
          fullranges.push_back (range_full_t{
              "positive normal", float_ranges_t<F>::plus_normal_min,
              float_ranges_t<F>::plus_normal_max });
          fullranges.push_back (range_full_t{
              "negative normal", float_ranges_t<F>::neg_normal_min,
              float_ranges_t<F>::neg_normal_max });
        }
      else if (ptrange == "subnormal")
        {
          fullranges.push_back (range_full_t{
              "positive subnormal", float_ranges_t<F>::plus_subnormal_min,
              float_ranges_t<F>::plus_subnormal_max });
          fullranges.push_back (range_full_t{
              "negative subnormal", float_ranges_t<F>::neg_subnormal_min,
              float_ranges_t<F>::neg_subnormal_max });
        }
    }
  return fullranges;
}

template <typename F>
const range_random_list_t<F>
get_as_range_random_list (const ranges_types_t<F> &ranges)
{
  return std::get<range_random_list_t<F> > (ranges);
}

template <typename F>
const range_full_list_t
get_as_range_full_list (const ranges_types_t<F> &ranges)
{
  return std::get<range_full_list_t> (ranges);
}

template <typename F>
static void
run_univariate (const std::string &function, bool coremath,
                const round_set &round_modes, fail_mode_t failmode,
                const boost::property_tree::ptree &jsontree, bool verbose)
{
  auto ranges = parse_ranges<F> (jsontree, verbose);
  auto func = get_univariate<F> (function, coremath).value ();
  if (!func.first)
    error ("libc does not provide {}", function);

  switch (ranges.index ())
    {
    case range_type_mode_t::RANGE_RANDOM:
      {
        random_univariate_t<F> sample{ func.first, func.second };
        check_random_variate (function, sample,
                              get_as_range_random_list<F> (ranges),
                              round_modes, failmode);
      }
      break;
    case range_type_mode_t::RANGE_FULL:
      {
        full_univariate_t<F> sample{ func.first, func.second };
        check_full_variate (function, sample,
                            get_as_range_full_list<F> (ranges), round_modes,
                            failmode);
      }
      break;
    }
}

template <typename F>
static void
run_bivariate (const std::string &function, bool coremath,
               const round_set &round_modes, fail_mode_t failmode,
               const boost::property_tree::ptree &jsontree, bool verbose)
{
  auto ranges = parse_ranges<F> (jsontree, verbose);
  auto func = get_bivariate<F> (function, coremath).value ();
  if (!func.first)
    error ("libc does not provide {}", function);

  switch (ranges.index ())
    {
    case range_type_mode_t::RANGE_RANDOM:
      {
        random_bivariate_t<F> sample{ func.first, func.second };
        check_random_variate (function, sample,
                              get_as_range_random_list<F> (ranges),
                              round_modes, failmode);
      }
      break;
    case range_type_mode_t::RANGE_FULL:
      // TODO
      break;
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
      = fail_mode_from_options (vm["fail"].as<std::string> ());

  std::ifstream file (vm["desc"].as<std::string> ());

  boost::property_tree::ptree jsontree;
  boost::property_tree::json_parser::read_json (file, jsontree);

  std::string function = jsontree.get<std::string> ("function");

  auto type = get_func_type (function);
  if (!type)
    error ("invalid function: {}", function);

  init_random_state ();

  switch (type.value ())
    {
    case refimpls::func_type_t::binary32_univariate:
      run_univariate<float> (function, coremath, round_modes, failmode,
                             jsontree, verbose);
      break;
    case refimpls::func_type_t::binary32_bivariate:
      run_bivariate<float> (function, coremath, round_modes, failmode,
                            jsontree, verbose);
      break;
    case refimpls::func_type_t::binary64_univariate:
      run_univariate<double> (function, coremath, round_modes, failmode,
                              jsontree, verbose);
      break;
    case refimpls::func_type_t::binary64_bivariate:
      run_bivariate<double> (function, coremath, round_modes, failmode,
                             jsontree, verbose);
      break;
    default:
      error ("function {} not implemented", function);
    }
}
