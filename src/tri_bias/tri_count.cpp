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

int CountTriNucleotide(const string& content, const string& codon) {
	string::size_type pos = 0;
	int count = 0;
	for (;;) {
	    pos = content.find(codon, pos);
	    if (pos == string::npos) {
	        break;
	    }
	    ++count;
	    ++pos;
	}
	return count;
}

vector<int> CountTriNucleotides(const string& content) {
	vector<future<int>> futures;
	for (auto codon : codons) {
		futures.push_back(async(CountTriNucleotide, content, codon));
	}
	vector<int> counts;
	for (auto& f : futures) {
		counts.push_back(f.get());
	}
	return counts;
}
