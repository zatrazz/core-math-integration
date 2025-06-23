#include <vector>
#include <algorithm>
#include <random>
#include <climits>
#include <iostream>
#include <ranges>
#include <string_view>

#include <boost/program_options.hpp>

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

static void
gen_range_binary32 (const std::string &name,
		    const std::string &start,
		    const std::string &end,
		    int count)
{
  float fstart = std::stof (start);
  float fend = std::stof (end);

  std::printf ("## args: float\n");
  std::printf ("## ret: float\n");
  std::printf ("## includes: math.h\n");
  std::printf ("## name: workload-%s\n", name.c_str());
  std::printf ("# Random inputs in [%.2f,%.2f]\n", fstart, fend);

  std::mt19937_64 e = seed_rng<std::mt19937_64>();
  std::uniform_real_distribution<float> d(fstart, fend);
  for (int i = 0; i < count; i++)
    std::printf ("%a\n", d(e));
}

static void
gen_range_binary64 (const std::string &name,
		    const std::string &start,
		    const std::string &end,
		    int count)
{
  double fstart = std::stod (start);
  double fend = std::stod (end);

  std::printf ("## args: double\n");
  std::printf ("## ret: double\n");
  std::printf ("## includes: math.h\n");
  std::printf ("## name: workload-%s\n", name.c_str());
  std::printf ("# Random inputs in [%.2lf,%.2lf]\n", fstart, fend);

  std::mt19937_64 e = seed_rng<std::mt19937_64>();
  std::uniform_real_distribution<double> d(fstart, fend);
  for (int i = 0; i < count; i++)
    std::printf ("%a\n", d(e));
}

static std::vector<std::string>
splitWithRanges(const std::string& s, std::string_view delimiter)
{
  std::vector<std::string> tokens;
  for (const auto& subrange : s | std::views::split(delimiter))
    tokens.push_back(std::string(subrange.begin(), subrange.end()));
  return tokens;
}

int main (int argc, char *argv[])
{
  namespace po = boost::program_options;

  po::options_description desc{"options"};
  desc.add_options()
    ("help,h",  "help")
    ("type,t",  po::value<std::string>()->default_value("binary32"), "type to use")
    ("range,r", po::value<std::string>()->default_value("0.0 10.0"), "range to use")
    ("count,c", po::value<int>()->default_value(1000)) 
    ("name,n",  po::value<std::string>());

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
  std::vector<std::string> rangeFields =
    splitWithRanges (range, " ");
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
    gen_range_binary64 (name, rangeFields[0], rangeFields[1], count);
  else if (type == "binary64")
    gen_range_binary64 (name, rangeFields[0], rangeFields[1], count);

  return 0;
}
