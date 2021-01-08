#pragma once
#include <iostream>
#include <string>
#include <unordered_map>
#include <boost/algorithm/string.hpp>
#include <xtensor/xarray.hpp>
#include <xtensor/xrandom.hpp>

using namespace std;

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
			vector<double> tempVector;
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

				if (tail == "right") {
					tempVector.push_back(stod(elements[1]));
				} else {
					tempVector.push_back(stod(elements[2]));
				}
				++i;
			}

			if (tempVector.empty()) {
				throw runtime_error("GeneProbabilies: no data");
			}

			probabilities = xt::zeros<double>({tempVector.size()});
			random_probabilities = xt::zeros<double>({tempVector.size()});
			for (int j = 0; j < tempVector.size(); ++j) {
				probabilities[j] = tempVector[j];
				random_probabilities[j] = tempVector[j];
			}
			xt::random::shuffle(random_probabilities);
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
