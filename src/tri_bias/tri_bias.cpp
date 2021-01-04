/*
Compute tri-nucleotide bias for assemblies in data/sequences. 
*/
#include <fstream>
#include <iostream>
#include <string>
#include <chrono>
#include <numeric>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include "../assembly/assembly.cpp"
#include "../fasta/fasta.cpp"
#include "../dna/dna.cpp"
#include "../count/tri_count.cpp"
#include "csv.hpp"

using namespace std;
using namespace std::chrono;
using namespace csv;

int main() {
	string dataFolder = "../data/";
	string sequencesFolder = dataFolder + "sequences/";
	string assembliesPath = dataFolder + "assemblies.csv";

	Assemblies assemblies(assembliesPath);

	cerr << "Processing " << assemblies.Size() << " assemblies" << endl;

	const bool overlap = true;

	vector<string> headers{"assembly_accession"};
	for (auto codon : codons) {
		headers.push_back(codon);
	}

	ofstream of("../data/tri_bias.csv");
	auto writer = make_csv_writer(of);
	writer << headers;

	int i = 0;
	auto start = system_clock::now();
	for (auto accession : assemblies.GetIds()) {
		if (i == 0 || (i+1) % 100 == 0) {
			auto tp = system_clock::now();
			cerr << "Processing assembly " << i + 1 << " / " << assemblies.Size();
			cerr << " (elapsed: " << duration_cast<seconds>(tp - start).count() << " seconds)" << endl;
		}
		string genome_path = sequencesFolder + accession + "/" + accession + "_genomic.fna.gz";
		ifstream input_file(genome_path, ios_base::in | ios_base::binary);
		boost::iostreams::filtering_streambuf<boost::iostreams::input> inbuf;
		inbuf.push(boost::iostreams::gzip_decompressor());
		inbuf.push(input_file);
		istream instream(&inbuf);

		FastaParser parser(instream);

		vector<int> totalCounts(codons.size());
		FastaRecord record{};
		while (parser.Get(record)) {
			auto counts = CountTriNucleotides(record.content, overlap);
			ReverseComplement(record.content);
			auto counts_rc = CountTriNucleotides(record.content, overlap);
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
	cerr << "DONE (elapsed: " << duration_cast<seconds>(tp - start).count() << " seconds)" << endl;
}
