#include <string>
#include <vector>
#include <future>

using namespace std;

const vector<string> codons = {
	"AAA", "AAC", "AAG", "AAT",
	"ACA", "ACC", "ACG", "ACT",
	"AGA", "AGC", "AGG", "AGT",
	"ATA", "ATC", "ATG", "ATT",
	"CAA", "CAC", "CAG", "CAT", 
	"CCA", "CCC", "CCG", "CCT",
	"CGA", "CGC", "CGG", "CGT",
	"CTA", "CTC", "CTG", "CTT",
	"GAA", "GAC", "GAG", "GAT",
	"GCA", "GCC", "GCG", "GCT",
	"GGA", "GGC", "GGG", "GGT",
	"GTA", "GTC", "GTG", "GTT",
	"TAA", "TAC", "TAG", "TAT",
	"TCA", "TCC", "TCG", "TCT",
	"TGA", "TGC", "TGG", "TGT",
	"TTA", "TTC", "TTG", "TTT",
};

int CountKmer(const string& content, const string& kmer, const bool overlap) {
	int count = 0;
	string::size_type pos = 0;
	for (;;) {
	    pos = content.find(kmer, pos);
	    if (pos == string::npos) {
	        break;
	    }
	    ++count;
	    if (overlap) {
	    	++pos;
	    } else {
	    	pos += kmer.size();
	    }
	}
	return count;
}

vector<int> CountTriNucleotidesAsync(const string& content, const bool overlap) {
	vector<future<int>> futures;
	for (auto& codon : codons) {
		futures.push_back(async(CountKmer, content, codon, overlap));
	}
	vector<int> counts;
	for (auto& f : futures) {
		counts.push_back(f.get());
	}
	return counts;
}

vector<int> CountTriNucleotidesSequential(const string& content, const bool overlap) {
	vector<int> counts;
	for (auto& codon : codons) {
		counts.push_back(CountKmer(content, codon, overlap));
	}
	return counts;
}

vector<int> CountTriNucleotides(const string& content, const bool overlap, const bool use_async) {
	if (use_async) {
		return CountTriNucleotidesAsync(content, overlap);
	} else {
		return CountTriNucleotidesSequential(content, overlap);
	}
}
