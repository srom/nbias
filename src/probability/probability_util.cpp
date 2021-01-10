#pragma once
#include <iostream>
#include <string>
#include <cmath>
#include <limits>
#include <unordered_map>
#include <boost/algorithm/string.hpp>
#include <xtensor/xarray.hpp>
#include <xtensor/xmath.hpp>

using namespace std;

xt::xarray<double> compute_probability_from_distance(const vector<double>& distances, const bool left_tail) {
	xt::xarray<double> probabilities = xt::zeros<double>({distances.size()});
	for (int i = 0; i < distances.size(); ++i) {
		auto d = distances[i];
		if (d > (double) 0) {
			probabilities[i] = log(distances[i]);
		} else {
			// Proxy for -infinity:
			probabilities[i] = (double) -10;
		}
	}

	double min_val = xt::amin(probabilities)[0];
	double max_val = xt::amax(probabilities)[0];

	double norm_min = floor(min_val);
	double norm_max = (double) 0;

	while (min_val <= norm_min) {
		norm_min -= 0.0001;
	}
	while (max_val >= norm_max) {
		norm_max += 0.0001;
	}

	if (left_tail) {
		return xt::eval((double) 1 - ((probabilities - norm_min) / (norm_max - norm_min)));
	} else {
		return xt::eval((probabilities - norm_min) / (norm_max - norm_min));
	}	
}

class GeneProbabilies {
	vector<string> keys;
	unordered_map<string, size_t> lookup;
	xt::xarray<double> probabilities;
	xt::xarray<double> random_probabilities;

	public:
		GeneProbabilies(istream& is, const string& tail) :keys{}, lookup{} {
			if (tail != "left" && tail != "right") {
				throw runtime_error(
					"GeneProbabilies: tail argument must be one of left, right "
					"but got: " + tail
				);
			}
			string line;
			int i = 0;
			vector<double> distances;
			while (getline(is, line)) {
				if (i == 0 || line.empty()) {
					++i;
					continue;
				}
				vector<string> elements;
				boost::split(elements, line, boost::is_any_of(","));
				if (elements.size() != 3) {
					throw runtime_error(
						"GeneProbabilies: 3 columns expected "
						"but got: " + to_string(elements.size())
					);
				}
				string key = elements[0];
				keys.push_back(key);
				lookup[key] = i-1;

				double distance = stod(elements[1]);
				distances.push_back(distance);
				++i;
			}

			if (distances.empty()) {
				throw runtime_error("GeneProbabilies: no data");
			}

			probabilities = compute_probability_from_distance(distances, tail == "left");

			double mean = xt::mean(probabilities)();
			random_probabilities = xt::zeros<double>({probabilities.size()});
			for (int j = 0; j < probabilities.size(); ++j) {
				random_probabilities[j] = mean;
			}
		}

		double operator[](const string& id) {
			return probabilities[lookup[id]];
		}

		double Get(const string& id, const bool random = false) {
			if (!random) {
				return probabilities[lookup[id]];
			} else {
				return random_probabilities[lookup[id]];
			}
		}

		size_t IndexOf(const string& id) {
			return lookup[id];
		}

		vector<string>& Keys() {
			return keys;
		}

		xt::xarray<double>& Values() {
			return probabilities;
		}

		xt::xarray<double>& RandomValues() {
			return random_probabilities;
		}

		size_t size() {
			return lookup.size();
		}
};

xt::xarray<double> make_uniform_prior(const size_t n_elements) {
	return xt::ones<double>({n_elements}) / n_elements;
}

xt::xarray<double> make_uniform_prior(const size_t n_rows, const size_t n_cols) {
	return xt::ones<double>({n_rows, n_cols}) / n_cols;
}
