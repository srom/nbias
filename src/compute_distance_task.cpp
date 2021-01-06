#include <fstream>
#include <iostream>
#include <string>
#include <chrono>
#include <future>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include "csv.hpp"
#include "assembly/assembly.cpp"
#include "distance/distance.cpp"
#include "distribution/distribution.cpp"

using namespace std;
using namespace std::chrono;
using namespace csv;

bool task_compute_distance(const int task_nb, const vector<string>& assembly_ids) {
	auto start = system_clock::now();

	cerr << "Thread " << task_nb << " started." << endl;

	string dataFolder = "../data/";
	string sequencesFolder = dataFolder + "sequences/";

	string assembly_dist_path = dataFolder + "tri_nucleotide_dist_genome_wide_with_overlap.csv";
	ifstream assembly_dist_f(assembly_dist_path);
	Distributions assembly_distributions(assembly_dist_f);

	int i = 0;
	for (auto& accession : assembly_ids) {
		if (i == 0 || (i+1) % 100 == 0) {
			auto tp = system_clock::now();
			cerr << "Thread " << task_nb << ": Processing assembly " << i + 1 << " / " << assembly_ids.size();
			cerr << " (elapsed: " << duration_cast<seconds>(tp - start).count() << " seconds)" << endl;
		}
		++i;

		string cds_dist_path = (
			sequencesFolder + accession + "/" + 
			accession + "_tri_nucleotide_dist_with_overlap.csv.gz"
		);
		string outputPath = (
			sequencesFolder + accession + "/" + 
			accession + "tri_nucleotide_distance_to_mean.csv"
		);

		ifstream cds_dist_f(cds_dist_path, ios_base::in | ios_base::binary);
		boost::iostreams::filtering_streambuf<boost::iostreams::input> inbuf;
		inbuf.push(boost::iostreams::gzip_decompressor());
		inbuf.push(cds_dist_f);
		istream instream(&inbuf);

		Distributions cds_distributions(instream);

		auto mean_dist = duplicate_rows(
			assembly_distributions[accession], 
			cds_distributions.size()
		);

		vector<string> protein_ids = cds_distributions.Keys();
		xt::xarray<double>& distributions = cds_distributions.Values();

		xt::xarray<double> distances = jensen_shannon_distance(mean_dist, distributions, 1);
		
		ofstream of(outputPath);
		auto writer = make_csv_writer(of);

		vector<string> headers{"protein_id", "distance", "probability"};
		writer << headers;

		for (int n = 0; n < distances.size(); ++n) {
			string protein_id = protein_ids[n];
			double distance = distances[n];
			double probability = 1.0 - distance;
			vector<string> row{protein_id, to_string(distance), to_string(probability)};
			writer << row;
		}
	}

	auto tp = system_clock::now();
	cerr << "Thread " << task_nb << ": DONE";
	cerr << " (elapsed: " << duration_cast<seconds>(tp - start).count() << " seconds)" << endl;
	return true;
}

void compute_distance(const int n_threads) {
	string dataFolder = "../data/";
	string assembliesPath = dataFolder + "assemblies.csv";

	Assemblies assemblies(assembliesPath);
	auto assembly_ids = assemblies.GetIds();

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
		futures.push_back(async(task_compute_distance, i+1, ids));
	}
	for (auto& f : futures) {
		if(!f.get()) {
			throw runtime_error("Unexpected error while processing future's output");
		}
	}
	cerr << "DONE" << endl;
}
