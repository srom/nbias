#include <fstream>
#include <iostream>
#include <filesystem>
#include <string>
#include <random>
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

unordered_map<string, vector<double>> MakeInterventions(
	FastaParser parser, 
	const xt::xarray<double>& genome_wide_dist,
	const string& intervention, 
	const int n_samples,
	mt19937 rng
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
				seq = SwapSynonymousCodons(sequence, rng);
			} else if (intervention == "shuffle") {
				seq = ShuffleCodons(sequence, rng);
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
				cerr << "MakeInterventions: Nan or Inf encountered for protein id: " << protein_id << endl;
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
	auto start_clock = system_clock::now();

	cerr << "Thread " << task_nb << " started." << endl;

	string dataFolder = "../data/";
	string sequencesFolder = dataFolder + "sequences/";
	auto n_assemblies = assembly_ids.size();

	// Load genome-wide distributions
	string assembly_dist_path = (
		dataFolder + "tri_nucleotide_dist_genome_wide_with_overlap.csv"
	);
	ifstream assembly_dist_f(assembly_dist_path);
	Distributions assembly_distributions(assembly_dist_f);

	// Load Protein domains metadata
	string metadata_path = dataFolder + query + "_master.csv";
	auto metadata = LoadDomainMetadata(metadata_path);

	// Random number generator
	random_device rd;
	mt19937 rng(rd());

	// Compute interventions on every assemblies
	int i = 0;
	for (auto& accession : assembly_ids) {
		auto elapsed = duration_cast<seconds>(system_clock::now() - start_clock).count();
		cerr << "Thread " << task_nb << ": ";
		cerr << "Processing assembly " << i + 1 << " / " << n_assemblies;
		cerr << " (elapsed: " << elapsed << " seconds)" << endl;
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
		unordered_map<string, vector<double>> protein_id_to_distances = MakeInterventions(
			parser,
			genome_wide_dist,
			intervention,
			n_samples,
			rng
		);

		// Compute mean probability per sample
		vector<double> means(n_samples);
		for (int s = 0; s < n_samples; ++s) {
			xt::xarray<double> probs = xt::zeros<double>({protein_id_to_distances.size()});
			int ix = 0;
			for (auto& it : protein_id_to_distances) {
				double distance = it.second[s];
				probs[ix] = compute_probability_from_distance(distance, tail == "left");
				++ix;
			}
			double mean = xt::mean(probs)();
			means[s] = mean;
		}

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
				double mean = means[s];

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

bool task_compute_domain_probabilities_per_phylum(
	const int task_nb,
	const string& query, 
	const string& tail,
	const string intervention,
	const vector<string>& phyla,
	const unordered_map<string, vector<string>>& assemblies_per_phylum
) {
	cerr << "Thread " << task_nb << " started." << endl;

	string dataFolder = "../data/intervention/" + intervention + "/";
	string sequencesFolder = dataFolder + "sequences/";
	string phylumFolder = dataFolder + "phylum/";
	auto n_phyla = phyla.size();

	// Create phylum directory if it does not exist.
	filesystem::create_directory(phylumFolder);

	int i = 0;
	for (auto& phylum : phyla) {
		auto& assembly_ids = assemblies_per_phylum.at(phylum);
		auto n_assemblies = assembly_ids.size();

		cerr << "Thread " << task_nb << ": ";
		cerr << "Processing phylum " << i + 1 << " / " << n_phyla;
		cerr << ": " << phylum << " (" << n_assemblies << " assemblies)";
		cerr << endl;
		++i;

		set<ProteinDomain> protein_domains;
		unordered_map<ProteinDomain, vector<DomainProbability>> protein_domain_probs;
		for (auto& accession : assembly_ids) {
			string path = (
				sequencesFolder + accession + "/" + 
				accession + "_" + query + "_probability_" + tail + ".csv"
			);
			vector<DomainProbability> domains = LoadDomainProbabilities(path);
			for (auto& domain_prob : domains) {
				auto& domain = domain_prob.domain;
				protein_domains.insert(domain);

				if (protein_domain_probs.find(domain) == protein_domain_probs.end()) {
					protein_domain_probs[domain] = vector<DomainProbability>{domain_prob};
				} else {
					protein_domain_probs[domain].push_back(domain_prob);
				}
			}
		}

		string phylum_lower = phylum;
		transform(phylum_lower.begin(), phylum_lower.end(), phylum_lower.begin(), ::tolower);
		transform(phylum_lower.begin(), phylum_lower.end(), phylum_lower.begin(), [](char ch) {
		    return ch == ' ' ? '_' : ch;
		});

		string phylumDir = phylumFolder + phylum_lower + "/";

		filesystem::create_directory(phylumDir);

		string phylum_domain_prob_out_path = (
			phylumDir + 
			phylum_lower + "_" + query + "_probability_" + tail + ".csv"
		);
		ofstream of(phylum_domain_prob_out_path);
		auto writer = make_csv_writer(of);
		writer << DomainProbability::RecordHeader();

		vector<DomainProbability> records;
		for (auto& domain : protein_domains) {
			auto& domain_probs = protein_domain_probs[domain];
			auto n_probs = domain_probs.size();
			xt::xarray<double> log_probs = xt::zeros<double>({n_probs});
			xt::xarray<double> log_probs_random = xt::zeros<double>({n_probs});
			for (int ix = 0; ix < n_probs; ++ix) {
				log_probs[ix] = domain_probs[ix].log_probability;
				log_probs_random[ix] = domain_probs[ix].log_probability_random;
			}

			double log_prob = product_rule_log(log_probs);
			double log_prob_random = product_rule_log(log_probs_random);

			try {
				DomainProbability record(
					domain, 
					log_prob, 
					log_prob_random,
					n_probs
				);
				records.push_back(record);
			}
			catch (exception& e) {
				cerr << "Thread " << task_nb << " | Phylum: " << phylum << " | ";
				cerr << "Exception: " << e.what() << endl;
				throw;
			}
		}

		sort(records.begin(), records.end(), greater<DomainProbability>()); 

		for (auto& record : records) {
			writer << record.Record();
		}
	}
	cerr << "Thread " << task_nb << ": DONE" << endl;
	return true;
}

void compute_intervention(
	const string query, 
	const string tail, 
	const string intervention,
	const int n_samples,
	const int n_threads
) {
	auto start_clock = system_clock::now();

	string dataFolder = "../data/";

	bool complete_genome_only = true;
	Assemblies assemblies(dataFolder + "intervention/" + "sample_assemblies.csv");
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

	// 
	// 2) Compute probability of domains for each phylum.
	//
	cerr << "Processing of domain probabilities per phylum" << endl;

	unordered_map<string, vector<string>> assemblies_per_phylum;
	for (auto& assembly_id : assembly_ids) {
		Assembly& assembly = assemblies.Get(assembly_id);
		string phylum = assembly.phylum;
		if (phylum.empty()) {
			continue;
		}
		if (assemblies_per_phylum.find(phylum) == assemblies_per_phylum.end()) {
			assemblies_per_phylum[phylum] = vector<string>{assembly_id};
		} else {
			assemblies_per_phylum[phylum].push_back(assembly_id);
		}
	}

	set<string> phyla_set;
	for (auto& assembly_id : assembly_ids) {
		Assembly& assembly = assemblies.Get(assembly_id);
		string phylum = assembly.phylum;
		if (!phylum.empty()) {
			phyla_set.insert(phylum);
		}
	}

	vector<string> phyla;
	phyla.assign(phyla_set.begin(), phyla_set.end());
	auto n_phyla = phyla.size();

	cerr << "Processing " << n_phyla << " phyla" << endl;
	cerr << "Starting " << n_threads << " threads" << endl;

	n_per_thread = ceil((double) n_phyla / (double) n_threads);

	vector<future<bool>> futuresP;
	for (int i = 0; i < n_threads; ++i) {
		auto start = phyla.begin() + i * n_per_thread;
		auto end = phyla.end();
		int endInt = i * n_per_thread + n_per_thread;
		if (endInt < n_phyla) {
			end = phyla.begin() + endInt;
		}
		futuresP.push_back(async(
			task_compute_domain_probabilities_per_phylum, 
			i+1, 
			query, 
			tail, 
			intervention,
			vector<string>(start, end),
			assemblies_per_phylum
		));
	}
	for (auto& f : futuresP) {
		if(!f.get()) {
			throw runtime_error("Unexpected error while processing phylum output");
		}
	}
	auto elapsed = duration_cast<seconds>(system_clock::now() - start_clock).count();
	cerr << "Processing of domain probabilities per phylum is complete" << endl;
	cerr << "Elapsed: " << elapsed << " seconds" << endl;

	//
	// 3) Compute global probability of domains.
	//
	cerr << "Processing of domain probabilities globally" << endl;
	set<ProteinDomain> protein_domains;
	unordered_map<ProteinDomain, vector<DomainProbability>> protein_domain_probs;
	for (auto& phylum : phyla) {
		string phylum_lower = phylum;
		transform(phylum_lower.begin(), phylum_lower.end(), phylum_lower.begin(), ::tolower);
		transform(phylum_lower.begin(), phylum_lower.end(), phylum_lower.begin(), [](char ch) {
		    return ch == ' ' ? '_' : ch;
		});
		string phylumDir = (
			dataFolder + "intervention/" + intervention + "/phylum/" + phylum_lower + "/"
		);
		string phylum_domain_prob_path = (
			phylumDir + 
			phylum_lower + "_" + query + "_probability_" + tail + ".csv"
		);

		vector<DomainProbability> domains = LoadDomainProbabilities(phylum_domain_prob_path);
		for (auto& domain_prob : domains) {
			auto& domain = domain_prob.domain;
			protein_domains.insert(domain);

			if (protein_domain_probs.find(domain) == protein_domain_probs.end()) {
				protein_domain_probs[domain] = vector<DomainProbability>{domain_prob};
			} else {
				protein_domain_probs[domain].push_back(domain_prob);
			}
		}
	}

	const string protein_out_path = (
		dataFolder + "intervention/" + intervention + "/" + 
		query + "_probability_" + tail + ".csv"
	);
	ofstream output_file(protein_out_path);
	auto writer = make_csv_writer(output_file);
	writer << DomainProbability::RecordHeader();

	xt::xarray<double> uniform_log_prior = xt::eval(
		xt::log(make_uniform_prior(n_phyla))
	);

	vector<DomainProbability> records;
	for (auto& domain : protein_domains) {
		auto& domain_probs = protein_domain_probs[domain];
		auto n_probs = domain_probs.size();
		xt::xarray<double> log_probs = xt::zeros<double>({n_phyla});
		xt::xarray<double> log_probs_random = xt::zeros<double>({n_phyla});
		for (int ix = 0; ix < n_probs; ++ix) {
			log_probs[ix] = domain_probs[ix].log_probability;
			log_probs_random[ix] = domain_probs[ix].log_probability_random;
		}

		double log_prob = marginalization_log(
			uniform_log_prior, 
			log_probs
		);
		double log_prob_random = marginalization_log(
			uniform_log_prior, 
			log_probs_random
		);

		try {
			DomainProbability record(
				domain, 
				log_prob, 
				log_prob_random,
				n_probs
			);
			records.push_back(record);
		}
		catch (exception& e) {
			cerr << "Global computation | ";
			cerr << "Exception: " << e.what() << endl;
			throw;
		}
	}

	sort(records.begin(), records.end(), greater<DomainProbability>()); 

	for (auto& record : records) {
		writer << record.Record();
	}

	auto elapsed_final = duration_cast<seconds>(system_clock::now() - start_clock).count();
	cerr << "Elapsed: " << elapsed_final << " seconds" << endl;
	cerr << "DONE" << endl;
}
