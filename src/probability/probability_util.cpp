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

double compute_probability_from_distance(const double distance, const bool left_tail) {
	double prob;
	if (left_tail) {
		prob = (double) 1.0 - distance;
	} else {
		prob = distance;
	}
	if (isnan(prob)) {
		throw runtime_error("compute_probability_from_distance: NaN values encountered");
	}
	return prob;
}

xt::xarray<double> compute_probabilities_from_distance(const vector<double>& distances, const bool left_tail) {
	xt::xarray<double> probabilities = xt::zeros<double>({distances.size()});
	for (int i = 0; i < distances.size(); ++i) {
		probabilities[i] = compute_probability_from_distance(distances[i], left_tail);
	}
	return probabilities;
}

xt::xarray<double> compute_baseline_probabilities(const xt::xarray<double>& probabilities) {
	xt::xarray<double> baseline = xt::zeros<double>({probabilities.size()});
	double mean = xt::mean(probabilities)();
	for (int j = 0; j < probabilities.size(); ++j) {
		baseline[j] = mean;
	}
	return baseline;
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
				if (elements.size() != 2) {
					throw runtime_error(
						"GeneProbabilies: 2 columns expected "
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

			probabilities = compute_probabilities_from_distance(distances, tail == "left");	
			random_probabilities = compute_baseline_probabilities(probabilities);
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

double normalize_probability(const double& probability, const size_t& n_elements) {
	return exp(log(probability) / n_elements);
}
