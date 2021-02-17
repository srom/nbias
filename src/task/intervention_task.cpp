#include <fstream>
#include <iostream>
#include <filesystem>
#include <string>
#include <unordered_map>
#include <algorithm>
#include <chrono>
#include <future>
#include <utility>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <xtensor/xarray.hpp>
#include <xtensor/xmath.hpp>
#include "csv.hpp"
#include "../assembly/assembly.cpp"
#include "../probability/probability.cpp"
#include "../probability/probability_util.cpp"
#include "../domain/domain.cpp"

using namespace std;
using namespace std::chrono;
using namespace csv;

void compute_intervention(
	const string query, 
	const string tail, 
	const string intervention,
	const int n_threads
) {

}
