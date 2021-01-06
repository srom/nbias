#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <unordered_map>
#include <boost/algorithm/string.hpp>
#include <xtensor/xarray.hpp>
#include <xtensor/xview.hpp>

using namespace std;

class Distributions {
	vector<string> keys;
	unordered_map<string, size_t> lookup;
	xt::xarray<double> values;

	public:
		Distributions(istream& is) :keys{}, lookup{} {
			string line;
			int i = 0;
			vector<xt::xarray<double>> tempVector;
			while (getline(is, line)) {
				if (i == 0 || line.empty()) {
					++i;
					continue;
				}
				vector<string> elements;
				boost::split(elements, line, boost::is_any_of(","));
				if (elements.size() < 2) {
					throw runtime_error("Invalid distribution data");
				}
				string key = elements[0];
				keys.push_back(key);
				lookup[key] = i-1;

				xt::xarray<double> distribution = xt::zeros<double>({elements.size() - 1});
				for (int k = 0; k < elements.size() - 1; ++k) {
					distribution[k] = stod(elements[k+1]);
				}
				tempVector.push_back(distribution);
				++i;
			}
			if (tempVector.empty()) {
				return;
			}

			size_t innerDim = tempVector[0].size();
			values = xt::zeros<double>({tempVector.size(), innerDim});
			for (int j = 0; j < tempVector.size(); ++j) {
				for (int l = 0; l < innerDim; ++l) {
					values(j, l) = tempVector[j][l];
				}
			}
		}

		xt::xarray<double> operator[](const string& id) {
			return xt::view(values, lookup[id]);
		}

		vector<string>& Keys() {
			return keys;
		}

		xt::xarray<double>& Values() {
			return values;
		}

		size_t size() {
			return lookup.size();
		}
};
