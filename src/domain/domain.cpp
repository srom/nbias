#include <iostream>
#include <cmath>
#include <string>
#include <unordered_map>
#include <boost/algorithm/string.hpp>
#include <xtensor/xarray.hpp>
#include "csv.hpp"
#include "../probability/probability_util.cpp"

using namespace std;
using namespace csv;

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
				ProteinDomain domain(elements[3], elements[2]);

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
	} else if (evidence <= 1) {
		return "Substantial";
	} else if (evidence <= 2) {
		return "Strong";
	} else {
		return "Decisive";
	}
}

struct DomainProbability {
	ProteinDomain domain;
	double probability;
	double probability_random;
	double evidence;
	string evidence_strength;

	DomainProbability(
		const ProteinDomain& d, 
		const double& p,
		const double& pr
	) {
		domain = ProteinDomain(d);
		probability = p;
		probability_random = pr;
		evidence = log10(probability) - log10(probability_random);
		evidence_strength = EvidenceStrength(evidence);
	}

	vector<string> Record() {
		return vector<string>{
			domain.id,
			domain.query,
			to_string(probability),
			to_string(probability_random),
			to_string(evidence),
			evidence_strength,
		};
	}

	static vector<string> RecordHeader() {
		return vector<string>{
			"id",
			"query",
			"probability",
			"probability_random",
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
				row["probability"].get<double>(),
				row["probability_random"].get<double>()
			)
		);
	}
	return out;
}
