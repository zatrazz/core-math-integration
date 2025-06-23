#include <vector>
#include <algorithm>
#include <random>
#include <climits>
#include <iostream>
#include <ranges>
#include <string_view>
#include <format>

#include <boost/program_options.hpp>

namespace po = boost::program_options;

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

template <typename T>
struct ftypeinfo
{
  static std::string name();
  static T fromstring (const std::string &);
};

template <>
struct ftypeinfo<float>
{
  using type = float;
  static std::string name() { return "float"; }
  static float fromstring (const std::string &s) { return std::stof (s); }
};

template <>
struct ftypeinfo<double>
{
  using type = double;
  static std::string name() { return "double"; }
  static double fromstring (const std::string &s) { return std::stod (s); }
};

template <typename ftypeinfo>
static void gen_range_binary (const std::string &name,
			      const std::string &start,
			      const std::string &end,
			      int count)
{
  using ftype = typename ftypeinfo::type;

  ftype fstart = ftypeinfo::fromstring (start);
  ftype fend =   ftypeinfo::fromstring (end);

  println ("## args: {}", ftypeinfo::name());
  println ("## ret: {}", ftypeinfo::name());
  println ("## includes: math.h");
  println ("## name: workload-{}", name);
  println ("# Random inputs in [{:.2f},{:.2f}]", fstart, fend);

  std::mt19937_64 e = seed_rng<std::mt19937_64>();
  std::uniform_real_distribution<ftype> d(fstart, fend);
  for (int i = 0; i < count; i++) {
    ftype f = d(e);
    bool isnegative = f < 0.0;
    println ("{}0x{:a}", isnegative ? "-": "", std::fabs(f));
  }
}

template <typename ftypeinfo>
static void gen_range_binary_bivariate (const std::string &name,
					const std::string &start0,
					const std::string &end0,
					const std::string &start1,
					const std::string &end1,
					int count)
{
  using ftype = typename ftypeinfo::type;

  ftype fstart0 = ftypeinfo::fromstring (start0);
  ftype fend0 =   ftypeinfo::fromstring (end0);
  ftype fstart1 = ftypeinfo::fromstring (start1);
  ftype fend1 =   ftypeinfo::fromstring (end1);

  println ("## args: {0}:{0}", ftypeinfo::name());
  println ("## ret: {}", ftypeinfo::name());
  println ("## includes: math.h");
  println ("## name: workload-{}", name);
  println ("# Random inputs with in in [{:.2f},{:.2f}] and y in [{:.2f},{:.2f}]",
	   fstart0, fend0,
	   fstart1, fend1);

  std::mt19937_64 e = seed_rng<std::mt19937_64>();
  std::uniform_real_distribution<ftype> d0(fstart0, fend0);
  std::uniform_real_distribution<ftype> d1(fstart1, fend1);
  for (int i = 0; i < count; i++)
    {
      ftype f0 = d0(e);
      ftype f1 = d1(e);
      bool isnegative0 = f0 < 0.0;
      bool isnegative1 = f1 < 0.0;
      println ("{}0x{:a}, {}0x{:a}",
	       isnegative0 ? "-": "", std::fabs(f0),
	       isnegative1 ? "-": "", std::fabs(f1));
    }
}

static void
gen_range_binary32 (const std::string &name,
		    const std::string &start,
		    const std::string &end,
		    int count)
{
  return gen_range_binary<ftypeinfo<float>>(name, start, end, count);
}

static void
gen_range_binary32_bivariate (const std::string &name,
			      const std::string &start0,
			      const std::string &end0,
			      const std::string &start1,
			      const std::string &end1,
			      int count)
{
  return gen_range_binary_bivariate<ftypeinfo<float>>(name, start0, end0, start1, end1, count);
}

static void
gen_range_binary64 (const std::string &name,
		    const std::string &start,
		    const std::string &end,
		    int count)
{
  return gen_range_binary<ftypeinfo<double>>(name, start, end, count);
}

static void
gen_range_binary64_bivariate (const std::string &name,
			      const std::string &start0,
			      const std::string &end0,
			      const std::string &start1,
			      const std::string &end1,
			      int count)
{
  return gen_range_binary_bivariate<ftypeinfo<double>>(name, start0, end0, start1, end1, count);
}

static std::vector<std::string>
splitWithRanges(const std::string& s, std::string_view delimiter)
{
  std::vector<std::string> tokens;
  for (const auto& subrange : s | std::views::split(delimiter))
    tokens.push_back(std::string(subrange.begin(), subrange.end()));
  return tokens;
}

static int
handle_univariate (const po::variables_map &vm,
		   const std::string &type,
		   int count,
		   std::string &range)
{
  std::vector<std::string> rangeFields =
    splitWithRanges (range, ":");
  if (rangeFields.size() != 2)
    {
      std::cout << "error: invalid range (" << range << ")\n";
      return 1;
    }
  
  std::string name;
  if (vm.count("name"))
    name = vm["name"].as<std::string>();
  else
    name = "random";

  if (type == "binary32")
    gen_range_binary32 (name, rangeFields[0], rangeFields[1], count);
  else if (type == "binary64")
    gen_range_binary64 (name, rangeFields[0], rangeFields[1], count);

  return 0;
}

static int
handle_bivariate (const po::variables_map &vm,
		  const std::string &type,
		  int count,
		  std::string &range)
{
  std::vector<std::string> ranges = splitWithRanges(range, " ");
  if (ranges.size() != 2)
    {
      std::cout << "error: invalid range (" << range << ")\n";
      return 1;
    }

  std::vector<std::string> range0 = splitWithRanges (ranges[0], ":");
  std::vector<std::string> range1 = splitWithRanges (ranges[1], ":");
  if (range0.size() != 2 || range1.size() != 2)
    {
      std::cout << "error: invalid range (" << range << ")\n";
      return 1;
    }

  std::string name;
  if (vm.count("name"))
    name = vm["name"].as<std::string>();
  else
    name = "random";

  if (type == "binary32")
    gen_range_binary32_bivariate (name,
				  range0[0],
				  range0[1],
				  range1[0],
				  range1[1],
				  count);
  else if (type == "binary64");
    gen_range_binary64_bivariate (name,
				  range0[0],
				  range0[1],
				  range1[0],
				  range1[1],
				  count);

  return 0;
}


int main (int argc, char *argv[])
{
  po::options_description desc{"options"};
  desc.add_options()
    ("help,h",      "help")
    ("type,t",      po::value<std::string>()->default_value("binary32"), "type to use")
    ("range,r",     po::value<std::string>()->default_value("0.0:10.0"), "range to use")
    ("count,c",     po::value<int>()->default_value(1000)) 
    ("name,n",      po::value<std::string>())
    ("bivariate,b", po::bool_switch()->default_value(false)); 

  po::variables_map vm;
  po::store (po::parse_command_line (argc, argv, desc), vm);
  po::notify (vm);

  if (vm.count("help"))
    {
      std::cout << desc << "\n";
      return 1;
    }

  std::string type  = vm["type"].as<std::string>();
  int count         = vm["count"].as<int>();
  std::string range = vm["range"].as<std::string>();
  bool bivrariate   = vm["bivariate"].as<bool>();

  if (!bivrariate)
    return handle_univariate (vm, type, count, range);
  else
    return handle_bivariate (vm, type, count, range);
}
