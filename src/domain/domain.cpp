#include <iostream>
#include <cmath>
#include <string>
#include <sstream>
#include <unordered_map>
#include <boost/algorithm/string.hpp>
#include <xtensor/xarray.hpp>
#include "csv.hpp"
#include "../probability/probability_util.cpp"

using namespace std;
using namespace csv;

const double LOG_10 = log((double) 10);

struct ProteinDomain {
	string id;
	string query;

	ProteinDomain() :id{""}, query{""} {}

	ProteinDomain(
		const string& domain_id, 
		const string& domain_query
	) {
		id = domain_id;
		query = domain_query;
	}

	ProteinDomain(const ProteinDomain& d) {
		id = d.id;
		query = d.query;
	}
};

bool operator==(const ProteinDomain& a, const ProteinDomain& b) {
	return a.id == b.id;
}

bool operator!=(const ProteinDomain& a, const ProteinDomain& b) {
	return !(a == b);
}

bool operator<(const ProteinDomain& a, const ProteinDomain& b) {
	return a.id < b.id;
}

bool operator<=(const ProteinDomain& a, const ProteinDomain& b) {
	return (a < b || a == b) && !(b < a);
}

bool operator>(const ProteinDomain& a, const ProteinDomain& b) {
	return b < a;
}

bool operator>=(const ProteinDomain& a, const ProteinDomain& b) {
	return (a > b || a == b) && !(a < b);
}

namespace std {
	template <>
  	struct hash<ProteinDomain> {
		size_t operator()(const ProteinDomain& pd) const {
		  	return hash<string>()(pd.id);
		}
  };
}

class ProteinDomains {
	vector<ProteinDomain> keys;
	unordered_map<ProteinDomain, vector<string>> protein_id_map;

	public:
		ProteinDomains(istream& is) :keys{}, protein_id_map{} {
			string line;
			int i = 0;
			while (getline(is, line)) {
				if (i == 0 || line.empty()) {
					++i;
					continue;
				}
				++i;

				vector<string> elements;
				boost::split(elements, line, boost::is_any_of(","));
				if (elements.size() < 4) {
					throw runtime_error(
						"ProteinDomains: at least 4 columns expected "
						"but got: " + to_string(elements.size())
					);
				}
				string protein_id = elements[1];
				string query = elements[2];
				string full_id = elements[3];

				vector<string> id_parts;
				boost::split(id_parts, full_id, boost::is_any_of("."));
				string id = id_parts[0];
				
				ProteinDomain domain(id, query);

				if (protein_id_map.find(domain) == protein_id_map.end()) {
					protein_id_map[domain] = vector<string>{protein_id};
					keys.push_back(domain);
				} else {
					protein_id_map[domain].push_back(protein_id);
				}
			}
		}

		vector<ProteinDomain> Keys() {
			return keys;
		}

		vector<string> ProteinIds(const ProteinDomain& key) {
			return protein_id_map[key];
		}

		size_t size() {
			return keys.size();
		}

		xt::xarray<double> Probabilities(
			const ProteinDomain& key,
			GeneProbabilies& probs, 
			const bool random = false
		) {
			auto protein_ids = ProteinIds(key);
			xt::xarray<double> probabilities = xt::zeros<double>({protein_ids.size()});
			int k = 0;
			for (const string protein_id : protein_ids) {
				probabilities[k] = probs.Get(protein_id, random);
				++k;
			}
			return probabilities;
		}
};

string EvidenceStrength(double evidence) {
	if (evidence <= 0) {
		return "Negative";
	} else if (evidence <= 0.5) {
		return "Weak";
	} else if (evidence <= 1.0) {
		return "Substantial";
	} else if (evidence <= 1.5) {
		return "Strong";
	} else if (evidence <= 2.0) {
		return "Very Strong";
	} else {
		return "Decisive";
	}
}

template <typename T>
string to_string_with_precision(const T val, const int n = 8)
{
    ostringstream out;
    out.precision(n);
    out << fixed << val;
    return out.str();
}

struct DomainProbability {
	ProteinDomain domain;
	double log_probability;
	double log_probability_random;
	size_t n_elements;
	double evidence;
	string evidence_strength;

	DomainProbability(
		const ProteinDomain& d, 
		const double& log_p,
		const double& log_pr,
		const size_t n_el
	) {
		if (isinf(log_p) || isinf(log_pr)) {
			throw runtime_error("DomainProbability: -inf log probability encountered");
		} else if (isnan(log_p) || isnan(log_pr)) {
			throw runtime_error("DomainProbability: NaN log probability encountered");
		} else if (log_p > 0 || log_pr > 0) {
			throw runtime_error("DomainProbability: log probability > 0 encountered");
		} else if (n_el == 0) {
			throw runtime_error("DomainProbability: n_elements is zero");
		}

		domain = ProteinDomain(d);
		log_probability = log_p;
		log_probability_random = log_pr;
		n_elements = n_el;
		evidence = (log_probability - log_probability_random) / LOG_10;
		evidence_strength = EvidenceStrength(evidence);
	}

	vector<string> Record() {
		return vector<string>{
			domain.id,
			domain.query,
			to_string_with_precision(log_probability),
			to_string_with_precision(log_probability_random),
			to_string(n_elements),
			to_string_with_precision(evidence),
			evidence_strength,
		};
	}

	static vector<string> RecordHeader() {
		return vector<string>{
			"id",
			"query",
			"log_probability",
			"log_probability_random",
			"n_elements",
			"evidence",
			"evidence_strength",
		};
	}
};

bool operator==(const DomainProbability& a, const DomainProbability& b) {
	return a.evidence == b.evidence;
}

bool operator!=(const DomainProbability& a, const DomainProbability& b) {
	return !(a == b);
}

bool operator<(const DomainProbability& a, const DomainProbability& b) {
	return a.evidence < b.evidence;
}

bool operator<=(const DomainProbability& a, const DomainProbability& b) {
	return (a < b || a == b) && !(b < a);
}

bool operator>(const DomainProbability& a, const DomainProbability& b) {
	return b < a;
}

bool operator>=(const DomainProbability& a, const DomainProbability& b) {
	return (a > b || a == b) && !(a < b);
}

vector<DomainProbability> LoadDomainProbabilities(const string& path) {
	CSVReader reader(path);
	vector<DomainProbability> out;
	for (CSVRow& row: reader) {
		ProteinDomain domain(
			row["id"].get<string>(),
			row["query"].get<string>()
		);
		out.push_back(
			DomainProbability(
				domain,
				row["log_probability"].get<double>(),
				row["log_probability_random"].get<double>(),
				row["n_elements"].get<size_t>()
			)
		);
	}
	return out;
}
