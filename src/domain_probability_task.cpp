#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>
#include <chrono>
#include <future>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <xtensor/xarray.hpp>
#include <xtensor/xrandom.hpp>
#include "csv.hpp"
#include "assembly/assembly.cpp"
#include "probability/probability.cpp"
#include "probability/probability_util.cpp"
#include "domain/domain.cpp"

using namespace std;
using namespace std::chrono;
using namespace csv;

bool task_compute_domain_probabilities_per_assembly(
	const int task_nb,
	const string query, 
	const string tail,
	const vector<string>& assembly_ids
) {
	auto start = system_clock::now();

	cerr << "Thread " << task_nb << " started." << endl;

	string dataFolder = "../data/";
	string sequencesFolder = dataFolder + "sequences/";
	auto n_assemblies = assembly_ids.size();

	int i = 0;
	for (auto& accession : assembly_ids) {
		if (i == 0 || (i+1) % 100 == 0) {
			auto tp = system_clock::now();
			auto elapsed = duration_cast<seconds>(tp - start).count();
			cerr << "Thread " << task_nb << ": ";
			cerr << "Processing assembly " << i + 1 << " / " << n_assemblies;
			cerr << " (elapsed: " << elapsed << " seconds)" << endl;
		}
		++i;

		string gene_probs_path = (
			sequencesFolder + accession + "/" + 
			accession + "_tri_nucleotide_distance_to_mean.csv"
		);
		ifstream gene_probs_file(gene_probs_path); 
		GeneProbabilies gene_probs(gene_probs_file, tail);

		string protein_domains_path = (
			sequencesFolder + accession + "/" + 
			accession + "_" + query + ".csv.gz"
		);
		ifstream protein_domains_file(protein_domains_path);
		boost::iostreams::filtering_streambuf<boost::iostreams::input> inbuf;
		inbuf.push(boost::iostreams::gzip_decompressor());
		inbuf.push(protein_domains_file);
		istream instream(&inbuf);
		ProteinDomains domains(instream);

		string assembly_domain_prob_out_path = (
			sequencesFolder + accession + "/" + 
			accession + "_" + query + "_probability_" + tail + ".csv"
		);
		ofstream of(assembly_domain_prob_out_path);
		auto writer = make_csv_writer(of);
		writer << DomainProbability::RecordHeader();

		vector<DomainProbability> records;
		for (const ProteinDomain& domain : domains.Keys()) {
			xt::xarray<double> probabilities = domains.Probabilities(domain, gene_probs);
			xt::xarray<double> random_probs = domains.Probabilities(domain, gene_probs, true);

			double probability = product_rule(probabilities);
			double probability_random = product_rule(random_probs);

			DomainProbability record(domain, probability, probability_random);
			records.push_back(record);
		}

		sort(records.begin(), records.end(), greater<DomainProbability>()); 

		for (auto& record : records) {
			writer << record.Record();
		}
	}

	auto tp = system_clock::now();
	auto elapsed = duration_cast<seconds>(tp - start).count();
	cerr << "Thread " << task_nb << ": DONE";
	cerr << " (elapsed: " << elapsed << " seconds)" << endl;
	return true;
}

void compute_domain_probabilities(
	const string query, 
	const string tail, 
	const int n_threads,
	const int seed
) {
	xt::random::seed(seed);

	string dataFolder = "../data/";
	string assembliesPath = dataFolder + "assemblies.csv";

	Assemblies assemblies(assembliesPath);
	auto assembly_ids = assemblies.GetIds();

	cerr << "Processing of domain probabilities per assembly" << endl;
	cerr << "Processing " << assemblies.Size() << " assemblies" << endl;
	cerr << "Starting " << n_threads << " threads" << endl;

	auto n_per_thread = ceil((double) assembly_ids.size() / (double) n_threads);

	vector<future<bool>> futures;
	for (int i = 0; i < n_threads; ++i) {
		auto start = assembly_ids.begin() + i * n_per_thread;
		auto end = assembly_ids.end();
		int endInt = i * n_per_thread + n_per_thread;
		if (endInt < assembly_ids.size()) {
			end = assembly_ids.begin() + endInt;
		}
		auto ids = vector<string>(start, end);
		futures.push_back(async(
			task_compute_domain_probabilities_per_assembly, 
			i+1, 
			query, 
			tail, 
			ids
		));
	}
	for (auto& f : futures) {
		if(!f.get()) {
			throw runtime_error("Unexpected error while processing future's output");
		}
	}
	cerr << "Processing of domain probabilities per assembly is complete" << endl;
}
