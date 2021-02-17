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

bool task_intervention_per_assembly(
	const int task_nb,
	const string& query, 
	const string& tail,
	const string intervention,
	const int n_samples,
	const vector<string>& assembly_ids,
	const time_point start
) {
	cerr << "Thread " << task_nb << " started." << endl;

	string dataFolder = "../data/";
	string sequencesFolder = dataFolder + "sequences/";
	auto n_assemblies = assembly_ids.size();

	string assembly_dist_path = dataFolder + "tri_nucleotide_dist_genome_wide_with_overlap.csv";
	ifstream assembly_dist_f(assembly_dist_path);
	Distributions assembly_distributions(assembly_dist_f);

	string metadata_path = dataFolder + query + "_master.csv";
	auto metadata = LoadDomainMetadata(metadata_path);

	int i = 0;
	for (auto& accession : assembly_ids) {
		if (i == 0 || (i+1) % 100 == 0) {
			auto elapsed = duration_cast<seconds>(system_clock::now() - start).count();
			cerr << "Thread " << task_nb << ": ";
			cerr << "Processing assembly " << i + 1 << " / " << n_assemblies;
			cerr << " (elapsed: " << elapsed << " seconds)" << endl;
		}
		++i;

		xt::xarray<double>& genome_wide_dist = assembly_distributions[accession];

		cds_path = (
			sequencesFolder + 
			accession + "/" + 
			accession + "_cds_from_genomic.fna.gz"
		);
		ifstream input_file(cds_path, ios_base::in | ios_base::binary);
		boost::iostreams::filtering_streambuf<boost::iostreams::input> inbuf;
		inbuf.push(boost::iostreams::gzip_decompressor());
		inbuf.push(input_file);
		istream instream(&inbuf);

		FastaParser parser(instream);

		FastaRecord record{};
		while (parser.Get(record)) {
			auto& sequence = record.content;

			// Perform a number of random interventions on the CDS
			for (int s = 0; s < n_samples; ++s) {
				string seq = "";
				if (intervention == 'swap') {
					seq = SwapSynonymousCodons(sequence);
				} else {
					seq = ShuffleCodons(sequence);
				}

				// Count trinucleotide triplets
				bool overlap = true;
				bool use_async = false;
				auto counts = CountTriNucleotides(seq, overlap, use_async);
			}
		}
	}
}

void compute_intervention(
	const string query, 
	const string tail, 
	const string intervention,
	const int n_samples,
	const int n_threads
) {
	auto start = system_clock::now();

	string dataFolder = "../data/";

	bool complete_genome_only = true;
	Assemblies assemblies(dataFolder + "assemblies.csv", complete_genome_only);
	auto assembly_ids = assemblies.GetIds();
	auto n_per_thread = ceil((double) assembly_ids.size() / (double) n_threads);

	cerr << "Processing " << assemblies.Size() << " assemblies" << endl;

	// 
	// 1) Proceed with intervention and compute domain probability per assembly.
	//
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
			task_intervention_per_assembly, 
			i+1, 
			query, 
			tail, 
			intervention,
			n_samples,
			ids,
			start
		));
	}
	for (auto& f : futures) {
		if(!f.get()) {
			throw runtime_error("Unexpected error while processing assembly output");
		}
	}

	auto elapsed_final = duration_cast<seconds>(system_clock::now() - start).count();
	cerr << "Elapsed: " << elapsed_final << " seconds" << endl;
	cerr << "DONE" << endl;
}
