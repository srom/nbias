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
#include "assembly/assembly.cpp"
#include "probability/probability.cpp"
#include "probability/probability_util.cpp"
#include "domain/domain.cpp"

using namespace std;
using namespace std::chrono;
using namespace csv;

bool task_compute_domain_probabilities_per_assembly(
	const int task_nb,
	const string& query, 
	const string& tail,
	const vector<string>& assembly_ids
) {
	auto start = system_clock::now();

	cerr << "Thread " << task_nb << " started." << endl;

	string dataFolder = "../data/";
	string sequencesFolder = dataFolder + "sequences/";
	auto n_assemblies = assembly_ids.size();

	string metadata_path = dataFolder + query + "_master.csv";
	auto metadata = LoadDomainMetadata(metadata_path);

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

		try {
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
			for (ProteinDomain& domain : domains.Keys()) {
				if (metadata.find(domain.id) != metadata.end()) {
					auto& [domain_query, domain_description] = metadata[domain.id];
					domain.query = domain_query;
					domain.description = domain_description;
				}

				xt::xarray<double> probs = domains.Probabilities(domain, gene_probs);
				xt::xarray<double> probs_random = domains.Probabilities(domain, gene_probs, true);

				xt::xarray<double> log_probabilities = xt::eval(xt::log(probs));
				xt::xarray<double> log_probabilities_random = xt::eval(xt::log(probs_random));

				double log_prob = product_rule_log(log_probabilities);
				double log_prob_random = product_rule_log(log_probabilities_random);

				DomainProbability record(
					domain, 
					log_prob, 
					log_prob_random, 
					log_probabilities.size()
				);
				records.push_back(record);
			}

			sort(records.begin(), records.end(), greater<DomainProbability>()); 

			for (auto& record : records) {
				writer << record.Record();
			}
		}
		catch (exception& e) {
			cerr << "Thread " << task_nb << " | Assembly: " << accession << " | ";
			cerr << "Exception: " << e.what() << endl;
			throw;
		}
	}

	auto tp = system_clock::now();
	auto elapsed = duration_cast<seconds>(tp - start).count();
	cerr << "Thread " << task_nb << ": DONE";
	cerr << " (elapsed: " << elapsed << " seconds)" << endl;
	return true;
}

bool task_compute_domain_probabilities_per_phylum(
	const int task_nb,
	const string& query, 
	const string& tail,
	const vector<string>& phyla,
	const unordered_map<string, vector<string>>& assemblies_per_phylum
) {
	cerr << "Thread " << task_nb << " started." << endl;

	string dataFolder = "../data/";
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

void compute_domain_probabilities(
	const string query, 
	const string tail, 
	const int n_threads
) {
	auto start = system_clock::now();

	string dataFolder = "../data/";
	string assembliesPath = dataFolder + "assemblies.csv";

	Assemblies assemblies(assembliesPath);
	auto assembly_ids = assemblies.GetIds();
	auto n_per_thread = ceil((double) assembly_ids.size() / (double) n_threads);

	// 
	// 1) Compute probability of domains for each assembly individually.
	//
	cerr << "Processing of domain probabilities per assembly" << endl;
	cerr << "Processing " << assemblies.Size() << " assemblies" << endl;
	cerr << "Starting " << n_threads << " threads" << endl;

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
			throw runtime_error("Unexpected error while processing assembly output");
		}
	}

	auto tp = system_clock::now();
	auto elapsed = duration_cast<seconds>(tp - start).count();
	cerr << "Processing of domain probabilities per assembly is complete" << endl;
	cerr << "Elapsed: " << elapsed << " seconds" << endl;

	// 
	// 2) Compute probability of domains for each phylum 
	//    with at least 10 assemblies within it.
	//
	size_t min_n_phyla = 10;
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
		if (!phylum.empty() && assemblies_per_phylum[phylum].size() >= min_n_phyla) {
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
			vector<string>(start, end),
			assemblies_per_phylum
		));
	}
	for (auto& f : futuresP) {
		if(!f.get()) {
			throw runtime_error("Unexpected error while processing phylum output");
		}
	}

	tp = system_clock::now();
	elapsed = duration_cast<seconds>(tp - start).count();
	cerr << "Processing of domain probabilities per phylum is complete" << endl;
	cerr << "Elapsed: " << elapsed << " seconds" << endl;

	//
	// 3) Compute probability of domains per superkingdom.
	//
	cerr << "Processing of domain probabilities per superkingdom" << endl;
	vector<string> superkingdoms;
	unordered_map<string, vector<string>> phyla_per_superkingdom;
	for (auto& phylum : phyla) {
		auto assembly_id = assemblies_per_phylum[phylum][0];
		Assembly& assembly = assemblies.Get(assembly_id);
		const string superkingdom_raw = assembly.domain;
		if (superkingdom_raw.empty()) {
			continue;
		}
		string superkingdom = superkingdom_raw;
		transform(
			superkingdom.begin(), 
			superkingdom.end(), 
			superkingdom.begin(), 
			::tolower
		);
		superkingdom[0] = toupper(superkingdom[0]);

		if (phyla_per_superkingdom.find(superkingdom) == phyla_per_superkingdom.end()) {
			phyla_per_superkingdom[superkingdom] = vector<string>{phylum};
			superkingdoms.push_back(superkingdom);
		} else {
			phyla_per_superkingdom[superkingdom].push_back(phylum);
		}
	}

	string superkingdom_folder = dataFolder + "superkingdom/";
	filesystem::create_directory(superkingdom_folder);

	for (auto& superkingdom : superkingdoms) {
		string superkingdom_lower = superkingdom;
		transform(
			superkingdom_lower.begin(), 
			superkingdom_lower.end(), 
			superkingdom_lower.begin(), 
			::tolower
		);
		
		string superkingdom_inner_folder = (
			superkingdom_folder + "/" + superkingdom_lower + "/"
		);
		filesystem::create_directory(superkingdom_inner_folder);

		auto& superkingdom_phyla = phyla_per_superkingdom[superkingdom];
		auto n_superkingdom_phyla = superkingdom_phyla.size();

		set<ProteinDomain> protein_domains;
		unordered_map<ProteinDomain, vector<DomainProbability>> protein_domain_probs;
		for (auto& phylum : superkingdom_phyla) {
			string phylum_lower = phylum;
			transform(
				phylum_lower.begin(), 
				phylum_lower.end(), 
				phylum_lower.begin(), 
				::tolower
			);
			transform(
				phylum_lower.begin(), 
				phylum_lower.end(), 
				phylum_lower.begin(), 
				[](char ch) {
				    return ch == ' ' ? '_' : ch;
				}
			);
			string phylumDir = dataFolder + "phylum/" + phylum_lower + "/";
			string phylum_domain_prob_path = (
				phylumDir + 
				phylum_lower + "_" + query + "_probability_" + tail + ".csv"
			);
			vector<DomainProbability> domains = LoadDomainProbabilities(
				phylum_domain_prob_path
			);
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

		const string superkingdom_out_path = (
			superkingdom_inner_folder + 
			superkingdom_lower + "_" + query + "_probability_" + tail + ".csv"
		);
		ofstream superkingdom_of(superkingdom_out_path);
		auto writer = make_csv_writer(superkingdom_of);
		writer << DomainProbability::RecordHeader();

		xt::xarray<double> uniform_log_prior = xt::eval(
			xt::log(make_uniform_prior(n_superkingdom_phyla))
		);

		vector<DomainProbability> records;
		for (auto& domain : protein_domains) {
			auto& domain_probs = protein_domain_probs[domain];
			auto n_probs = domain_probs.size();
			xt::xarray<double> log_probs = xt::zeros<double>({n_superkingdom_phyla});
			xt::xarray<double> log_probs_random = xt::zeros<double>({n_superkingdom_phyla});
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
	}
	cerr << "Processing of domain probabilities per superkingdom is complete" << endl;

	//
	// 4) Compute global probability of domains.
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
		string phylumDir = dataFolder + "phylum/" + phylum_lower + "/";
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

	const string protein_out_path = dataFolder + query + "_probability_" + tail + ".csv";
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

	tp = system_clock::now();
	elapsed = duration_cast<seconds>(tp - start).count();
	cerr << "Processing of global domain probabilities is complete" << endl;
	cerr << "Elapsed: " << elapsed << " seconds" << endl;
	cerr << "DONE" << endl;
}
