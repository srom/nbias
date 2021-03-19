#pragma once
#include <fstream>
#include <string>
#include <json.hpp>

using namespace std;
using json = nlohmann::json;

struct DomainProbabilityContext {
	// Command line parameters
	const string kind;
	const string query;
	const string tail;
	const int n_threads;

	// Extra config
	bool complete_genome_only;
	string distance_to_mean_suffix;
	string assembly_output_folder;
	string phylum_output_folder;
	string superkingdom_output_folder;
	string overall_output_folder;

	DomainProbabilityContext(
		const string k,
		const string q,
		const string t,
		const int th,
		const bool cpo,
		const string dtm,
		const string assembly_folder,
		const string phylum_folder,
		const string superkingdom_folder,
		const string overall_output_folder
	) : kind(k), 
		query(q), 
		tail(t), 
		n_threads(th), 
		complete_genome_only(cpo),
		distance_to_mean_suffix(dtm),
		assembly_output_folder(assembly_folder),
		phylum_output_folder(phylum_folder),
		superkingdom_output_folder(superkingdom_folder),
		overall_output_folder(overall_output_folder)
	{

	}
};

DomainProbabilityContext parse_domain_probability_context(const string& path) {
	ifstream in(path);
	json j;
	in >> j;
	return DomainProbabilityContext{
		j["kind"].get<string>(),
		j["query"].get<string>(),
		j["tail"].get<string>(),
		j["n_threads"].get<int>(),
		j["complete_genome_only"].get<bool>(),
		j["distance_to_mean_suffix"].get<string>(),
		j["assembly_output_folder"].get<string>(),
		j["phylum_output_folder"].get<string>(),
		j["superkingdom_output_folder"].get<string>(),
		j["overall_output_folder"].get<string>()
	};
}

DomainProbabilityContext get_context_from_parameters(
	const string kind,
	const string query,
	const string tail,
	const int n_threads
) {
	string dtm_suffix;
	if (kind == "amino-acid") {
		dtm_suffix = "_amino_acid_distance_to_mean.csv";
	} else {
		dtm_suffix = "_tri_nucleotide_distance_to_mean.csv";
	}
	return DomainProbabilityContext{
		kind,
		query,
		tail,
		n_threads,
		false,
		dtm_suffix,
		"",
		"",
		"",
		""
	};
}
