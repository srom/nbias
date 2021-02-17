#include <fstream>
#include <iostream>
#include <filesystem>
#include <string>
#include <unordered_map>
#include <algorithm>
#include <chrono>
#include <future>
#include <utility>
#include <cmath>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <xtensor/xarray.hpp>
#include <xtensor/xmath.hpp>
#include "csv.hpp"
#include "../assembly/assembly.cpp"
#include "../distance/distance.cpp"
#include "../distribution/distribution.cpp"
#include "../probability/probability.cpp"
#include "../probability/probability_util.cpp"
#include "../domain/domain.cpp"
#include "../fasta/fasta.cpp"
#include "../count/tri_count.cpp"
#include "../codons/swap_codons.cpp"

using namespace std;
using namespace std::chrono;
using namespace csv;

unordered_map<string, vector<double>> MakeIntervention(
	FastaParser parser, 
	const xt::xarray<double>& genome_wide_dist,
	const string& intervention, 
	const int n_samples
) {
	unordered_map<string, vector<double>> protein_id_to_distances;
	FastaRecord cds_record{};
	while (parser.Get(cds_record)) {
		auto& sequence = cds_record.content;
		string protein_id = cds_record.id;

		// Perform a number of random interventions on the CDS
		vector<double> distances;
		for (int s = 0; s < n_samples; ++s) {
			string seq = "";
			if (intervention == "swap") {
				seq = SwapSynonymousCodons(sequence);
			} else if (intervention == "shuffle") {
				seq = ShuffleCodons(sequence);
			} else {
				throw runtime_error("Unknown intervention type: " + intervention);
			}

			// Count tri-nucleotides
			bool overlap = true;
			bool use_async = false;
			vector<int> counts = CountTriNucleotides(seq, overlap, use_async);

			// Compute distribution
			xt::xarray<double> cds_distribution = xt::zeros<double>({counts.size()});
			int sum = accumulate(counts.begin(), counts.end(), 0);
			for (int cix = 0; cix < counts.size(); ++cix) {
				cds_distribution[cix] = (double) counts[cix] / (double) sum;
			}

			// Compute distance
			double distance = jensen_shannon_distance(genome_wide_dist, cds_distribution);
			if (isinf(distance) || isnan(distance)) {
				cerr << "MakeIntervention: Nan or Inf encountered for protein id: " << protein_id << endl;
				continue;
			} else {
				distances.push_back(distance);
			}
		}
		if (distances.size() > 0) {
			protein_id_to_distances[protein_id] = distances;
		}
	}
	return protein_id_to_distances;
}

bool task_intervention_per_assembly(
	const int task_nb,
	const string& query, 
	const string& tail,
	const string intervention,
	const int n_samples,
	const vector<string>& assembly_ids
) {
	auto start = system_clock::now();

	cerr << "Thread " << task_nb << " started." << endl;

	string dataFolder = "../data/";
	string sequencesFolder = dataFolder + "sequences/";
	auto n_assemblies = assembly_ids.size();

	// Load genome-wide distributions
	string assembly_dist_path = dataFolder + "tri_nucleotide_dist_genome_wide_with_overlap.csv";
	ifstream assembly_dist_f(assembly_dist_path);
	Distributions assembly_distributions(assembly_dist_f);

	// Load Protein domains metadata
	string metadata_path = dataFolder + query + "_master.csv";
	auto metadata = LoadDomainMetadata(metadata_path);

	// Compute interventions on every assemblies
	int i = 0;
	for (auto& accession : assembly_ids) {
		if (i == 0 || (i+1) % 100 == 0) {
			auto elapsed = duration_cast<seconds>(system_clock::now() - start).count();
			cerr << "Thread " << task_nb << ": ";
			cerr << "Processing assembly " << i + 1 << " / " << n_assemblies;
			cerr << " (elapsed: " << elapsed << " seconds)" << endl;
		}
		++i;

		// Get genome-wide distribution
		xt::xarray<double> genome_wide_dist = assembly_distributions[accession];

		// Load protein domains to protein ids map
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

		// Instanciate CDS parser
		string cds_path = (
			sequencesFolder + 
			accession + "/" + 
			accession + "_cds_from_genomic.fna.gz"
		);
		ifstream input_file(cds_path, ios_base::in | ios_base::binary);
		boost::iostreams::filtering_streambuf<boost::iostreams::input> inbuf_cds;
		inbuf_cds.push(boost::iostreams::gzip_decompressor());
		inbuf_cds.push(input_file);
		istream instream_cds(&inbuf_cds);
		FastaParser parser(instream_cds);

		// Make interventions: randomly modify sequences a bunch of times and return results as a
		// map of protein ids to the list of distances from the modified sequence distributions 
		// the genome-wide distribution. 
		unordered_map<string, vector<double>> protein_id_to_distances = MakeIntervention(
			parser,
			genome_wide_dist,
			intervention,
			n_samples
		);

		// Compute domain probability
		vector<DomainProbability> records;
		for (ProteinDomain& domain : domains.Keys()) {
			// Set description from metadata rather than from annotations
			if (metadata.find(domain.id) != metadata.end()) {
				auto& [domain_query, domain_description] = metadata[domain.id];
				domain.query = domain_query;
				domain.description = domain_description;
			}

			// Get protein ids containing this domain
			auto protein_ids = domains.ProteinIds(domain);

			// Compute gene probabilities for each sample
			int n_records = -1;
			xt::xarray<double> log_probs = xt::zeros<double>({n_samples});
			xt::xarray<double> log_probs_baseline = xt::zeros<double>({n_samples});
			for (int s = 0; s < n_samples; ++s) {
				// Compute mean probability
				xt::xarray<double> probs = xt::zeros<double>({protein_id_to_distances.size()});
				int ix = 0;
				for (auto& it : protein_id_to_distances) {
					double distance = it.second[s];
					probs[ix] = compute_probability_from_distance(distance, tail == "left");
					++ix;
				}
				double mean = xt::mean(probs)();

				xt::xarray<double> probabilities = xt::zeros<double>({protein_ids.size()});
				xt::xarray<double> probabilities_baseline = xt::zeros<double>({protein_ids.size()});
				for (int k = 0; k < protein_ids.size(); ++k) {
					auto& protein_id = protein_ids[k];
					double distance = mean; // Default in case there is no data for protein id
					if (protein_id_to_distances.find(protein_id) != protein_id_to_distances.end()) {
						distance = protein_id_to_distances[protein_id][s];
					}
					probabilities[k] = compute_probability_from_distance(distance, tail == "left");
					probabilities_baseline[k] = mean;
				}

				xt::xarray<double> log_probabilities = xt::eval(xt::log(probabilities));
				xt::xarray<double> log_probabilities_baseline = xt::eval(xt::log(probabilities_baseline));

				double log_prob = product_rule_log(log_probabilities);
				double log_prob_baseline = product_rule_log(log_probabilities_baseline);

				if (n_records == -1) {
					n_records = log_probabilities.size();
				}
				log_probs[s] = log_prob;
				log_probs_baseline[s] = log_prob_baseline;
			}

			// Take the average of all samples
			DomainProbability record(
				domain, 
				xt::mean(log_probs)(), 
				xt::mean(log_probs_baseline)(), 
				n_records
			);
			records.push_back(record);
		}

		// Sort records from best to worse evidence
		sort(records.begin(), records.end(), greater<DomainProbability>()); 

		// Create assembly directory if it does not exist
		string assembly_domain_prob_out_folder = (
			dataFolder + "intervention/" +  intervention + "/sequences/" + accession + "/"
		);
		filesystem::create_directory(assembly_domain_prob_out_folder);

		// Prepare output writer
		string assembly_domain_prob_out_path = (
			assembly_domain_prob_out_folder +
			accession + "_" + query + "_probability_" + tail + ".csv"
		);
		ofstream of(assembly_domain_prob_out_path);
		auto writer = make_csv_writer(of);
		writer << DomainProbability::RecordHeader();

		// Write assembly outputs
		for (auto& record : records) {
			writer << record.Record();
		}
	}
	return true;
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
			ids
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
