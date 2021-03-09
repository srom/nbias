#include <fstream>
#include <iostream>
#include <string>
#include <chrono>
#include <numeric>
#include <cmath>
#include <cstdio>
#include <algorithm>
#include <boost/filesystem.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include "csv.hpp"
#include "../assembly/assembly.cpp"
#include "../fasta/fasta.cpp"
#include "../dna/dna.cpp"
#include "../count/tri_count.cpp"

using namespace std;
using namespace std::chrono;
using namespace csv;

string task_genome_wide_count(
	const int task_nb, 
	const string kind,
	const bool overlap, 
	const bool reverse_complement, 
	const vector<string>& assembly_ids
) {
	string dataFolder = "../data/";
	string sequencesFolder = dataFolder + "sequences/";

	auto p = boost::filesystem::temp_directory_path() / boost::filesystem::unique_path();
	string tempPath = p.native();

	cerr << "Thread " << task_nb << " started. Saving temp file to " << tempPath << endl;

	ofstream of(tempPath);
	auto writer = make_csv_writer(of);

	int i = 0;
	auto start = system_clock::now();
	bool use_async = true;
	for (auto& accession : assembly_ids) {
		if (i == 0 || (i+1) % 100 == 0) {
			auto tp = system_clock::now();
			cerr << "Thread " << task_nb << ": Processing assembly " << i + 1 << " / " << assembly_ids.size();
			cerr << " (elapsed: " << duration_cast<seconds>(tp - start).count() << " seconds)" << endl;
		}

		string genome_path;
		if (kind == "tri-nucleotide") {
			genome_path = sequencesFolder + accession + "/" + accession + "_genomic.fna.gz";
		} else {
			genome_path = sequencesFolder + accession + "/" + accession + "_protein.faa.gz";
		}

		ifstream input_file(genome_path, ios_base::in | ios_base::binary);
		boost::iostreams::filtering_streambuf<boost::iostreams::input> inbuf;
		inbuf.push(boost::iostreams::gzip_decompressor());
		inbuf.push(input_file);
		istream instream(&inbuf);

		FastaParser parser(instream);

		int size = 0;
		if (kind == "tri-nucleotide") {
			size = CODONS.size();
		} else {
			size = AMINO_ACIDS_THREE_LETTERS.size();
		}

		vector<int> totalCounts(size);
		FastaRecord record{};
		while (parser.Get(record)) {
			vector<int> counts(size);
			if (kind == "tri-nucleotide") {
				counts = CountTriNucleotides(record.content, overlap, use_async);
			} else {
				string content = record.content;
				transform(content.begin(), content.end(), content.begin(), ::toupper);
				counts = CountAminoAcids(content, false);
			}

			vector<int> counts_rc(size);
			if (reverse_complement) {
				ReverseComplement(record.content);
				counts_rc = CountTriNucleotides(record.content, overlap, use_async);
			}
			for (int k = 0; k < totalCounts.size(); ++k) {
				totalCounts[k] += counts[k] + counts_rc[k];
			}
		}
		int sum = accumulate(totalCounts.begin(), totalCounts.end(), 0);
		vector<string> row{accession};
		for (auto count : totalCounts) {
			row.push_back(to_string((double) count / (double) sum));
		}
		writer << row;
		++i;
	}
	auto tp = system_clock::now();
	cerr << "Thread " << task_nb << ": DONE";
	cerr << " (elapsed: " << duration_cast<seconds>(tp - start).count() << " seconds)" << endl;
	return tempPath;
}

void run_genome_wide_count(
	const string kind,
	const bool overlap,
	const bool reverse_complement, 
	const int n_threads
) {
	string dataFolder = "../data/";
	string assembliesPath = dataFolder + "assemblies.csv";

	Assemblies assemblies(assembliesPath);
	auto assembly_ids = assemblies.GetIds();

	cerr << "Processing " << assemblies.Size() << " assemblies" << endl;
	cerr << "Starting " << n_threads << " threads" << endl;

	auto n_per_thread = ceil((double) assembly_ids.size() / (double) n_threads);

	vector<future<string>> futurePaths;
	for (int i = 0; i < n_threads; ++i) {
		auto start = assembly_ids.begin() + i * n_per_thread;
		auto end = assembly_ids.end();
		int endInt = i * n_per_thread + n_per_thread;
		if (endInt < assembly_ids.size()) {
			end = assembly_ids.begin() + endInt;
		}
		auto ids = vector<string>(start, end);
		futurePaths.push_back(async(
			task_genome_wide_count, 
			i+1,  
			kind,
			overlap, 
			reverse_complement,
			ids
		));
	}

	vector<string> resultPaths;
	for (auto& f : futurePaths) {
		resultPaths.push_back(f.get());
	}

	string outputPath;

	if (kind == "tri-nucleotide") {
		outputPath = "../data/tri_nucleotide_dist_genome_wide";
	} else {
		outputPath = "../data/amino_acid_dist_genome_wide";
	}
	
	if (kind == "tri-nucleotide") {
		if (!reverse_complement) {
			outputPath += "_without_rc";
		}
		if (overlap) {
			outputPath += "_with_overlap";
		} else {
			outputPath += "_without_overlap";
		}
	}
	outputPath += ".csv";

	cerr << "Writing output file to " << outputPath << endl;

	ofstream writer(outputPath, ios_base::out);

	string header{"assembly_accession"};
	if (kind == "tri-nucleotide") {
		for (auto& codon : CODONS) {
			header += "," + codon;
		}
	} else {
		for (auto& aa_three_leters : AMINO_ACIDS_THREE_LETTERS) {
			header += "," + aa_three_leters;
		}
	}
	writer << header << endl;

	for (auto& path : resultPaths) {
		ifstream tempFile(path);
		if (!tempFile.good()) {
			throw runtime_error("Cannot open file " + path);
		}
		string line;
		while(getline(tempFile, line)) {
			writer << line << endl;
		}
		auto status = remove(path.c_str());
		if (status != 0) {
			cerr << "WARNING: Could not delete temp file " << path << endl;
		}
	}
	cerr << "DONE" << endl;
}

bool task_cds_count(
	const int task_nb, 
	const string kind,
	const bool overlap, 
	const bool reverse_complement, 
	const vector<string>& assembly_ids
) {
	string dataFolder = "../data/";
	string sequencesFolder = dataFolder + "sequences/";

	cerr << "Thread " << task_nb << " started." << endl;

	int i = 0;
	auto start = system_clock::now();
	for (auto& accession : assembly_ids) {
		if (i == 0 || (i+1) % 100 == 0) {
			auto tp = system_clock::now();
			cerr << "Thread " << task_nb << ": Processing assembly " << i + 1 << " / " << assembly_ids.size();
			cerr << " (elapsed: " << duration_cast<seconds>(tp - start).count() << " seconds)" << endl;
		}

		string cds_path;
		if (kind == "tri-nucleotide") {
			cds_path = (
				sequencesFolder + accession + "/" + 
				accession + "_cds_from_genomic.fna.gz"
			);
		} else {
			cds_path = (
				sequencesFolder + accession + "/" + 
				accession + "_protein.faa.gz"
			);
		}

		ifstream input_file(cds_path, ios_base::in | ios_base::binary);
		boost::iostreams::filtering_streambuf<boost::iostreams::input> inbuf;
		inbuf.push(boost::iostreams::gzip_decompressor());
		inbuf.push(input_file);
		istream instream(&inbuf);

		FastaParser parser(instream);

		string output_path;
		if (kind == "tri-nucleotide") {
			output_path = (
				sequencesFolder + accession + "/" + 
				accession + "_tri_nucleotide_dist"
			);
			if (!reverse_complement) {
				output_path += "_without_rc";
			}
			if (overlap) {
				output_path += "_with_overlap";
			} else {
				output_path += "_without_overlap";
			}
		} else {
			output_path = (
				sequencesFolder + accession + "/" + 
				accession + "_amino_acid_dist"
			);
		}
		output_path += ".csv.gz";

		ofstream of(output_path);
		boost::iostreams::filtering_streambuf<boost::iostreams::output> obuf;
		obuf.push(boost::iostreams::gzip_compressor());
		obuf.push(of);

		ostream ostream(&obuf);
		auto writer = make_csv_writer(ostream);

		vector<string> headers{"protein_id"};
		if (kind == "tri-nucleotide") {
			for (auto& codon : CODONS) {
				headers.push_back(codon);
			}
		} else {
			for (auto& aa : AMINO_ACIDS_THREE_LETTERS) {
				headers.push_back(aa);
			}
		}
		writer << headers;

		FastaRecord record{};
		while (parser.Get(record)) {
			vector<int> counts;
			if (kind == "tri-nucleotide") {
				counts = CountTriNucleotides(record.content, overlap, false);
			} else {
				string content = record.content;
				transform(content.begin(), content.end(), content.begin(), ::toupper);
				counts = CountAminoAcids(content, false);
			}

			if (reverse_complement) {
				ReverseComplement(record.content);
				auto counts_rc = CountTriNucleotides(record.content, overlap, false);
				for (int k = 0; k < counts.size(); ++k) {
					counts[k] += counts_rc[k];
				}
			}

			int sum = accumulate(counts.begin(), counts.end(), 0);
			vector<string> row{record.id};
			for (auto count : counts) {
				row.push_back(to_string((double) count / (double) sum));
			}
			writer << row;
		}
		++i;
	}
	auto tp = system_clock::now();
	cerr << "Thread " << task_nb << ": DONE";
	cerr << " (elapsed: " << duration_cast<seconds>(tp - start).count() << " seconds)" << endl;
	return true;
}

void run_cds_count(
	const string kind,
	const bool overlap, 
	const bool reverse_complement, 
	const int n_threads
) {
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
		futures.push_back(async(task_cds_count, i+1, kind, overlap, reverse_complement, ids));
	}
	for (auto& f : futures) {
		if(!f.get()) {
			throw runtime_error("Unexpected error while processing future's output");
		}
	}
	cerr << "DONE" << endl;
}
