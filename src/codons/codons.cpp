#include <string>
#include <vector>
#include <unordered_map>

using namespace std;

const vector<string> CODONS = {
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

const unordered_map<string, string> CODON_TO_AMINO_ACID = {
    {"TTT", "Phe"}, {"TTC", "Phe"}, {"TTA", "Leu"}, {"TTG", "Leu"},
    {"TCT", "Ser"}, {"TCC", "Ser"}, {"TCA", "Ser"}, {"TCG", "Ser"},
    {"TAT", "Tyr"}, {"TAC", "Tyr"}, {"TAA", "STP"}, {"TAG", "STP"},
    {"TGT", "Cys"}, {"TGC", "Cys"}, {"TGG", "Trp"}, {"TGA", "STP"},
    {"CTT", "Leu"}, {"CTC", "Leu"}, {"CTA", "Leu"}, {"CTG", "Leu"},
    {"CCT", "Pro"}, {"CCC", "Pro"}, {"CCA", "Pro"}, {"CCG", "Pro"},
    {"CAT", "His"}, {"CAC", "His"}, {"CAA", "Gln"}, {"CAG", "Gln"},
    {"CGT", "Arg"}, {"CGC", "Arg"}, {"CGA", "Arg"}, {"CGG", "Arg"},
    {"ATT", "Ile"}, {"ATC", "Ile"}, {"ATA", "Ile"}, {"ATG", "Met"},
    {"ACT", "Thr"}, {"ACC", "Thr"}, {"ACA", "Thr"}, {"ACG", "Thr"},
    {"AAT", "Asn"}, {"AAC", "Asn"}, {"AAA", "Lys"}, {"AAG", "Lys"},
    {"AGT", "Ser"}, {"AGC", "Ser"}, {"AGA", "Arg"}, {"AGG", "Arg"},
    {"GTT", "Val"}, {"GTC", "Val"}, {"GTA", "Val"}, {"GTG", "Val"},
    {"GCT", "Ala"}, {"GCC", "Ala"}, {"GCA", "Ala"}, {"GCG", "Ala"},
    {"GAT", "Asp"}, {"GAC", "Asp"}, {"GAA", "Glu"}, {"GAG", "Glu"},
    {"GGT", "Gly"}, {"GGC", "Gly"}, {"GGA", "Gly"}, {"GGG", "Gly"},
};

const unordered_map<string, vector<string>> AMINO_ACID_TO_CODONS = {
    {"Cys", {"TGT", "TGC"}},
    {"Asp", {"GAT", "GAC"}},
    {"Ser", {"TCT", "TCG", "TCA", "TCC", "AGC", "AGT"}},
    {"Gln", {"CAA", "CAG"}},
    {"Met", {"ATG"}},
    {"Asn", {"AAC", "AAT"}},
    {"Pro", {"CCT", "CCG", "CCA", "CCC"}},
    {"Lys", {"AAG", "AAA"}},
    {"STP", {"TAG", "TGA", "TAA"}},
    {"Thr", {"ACC", "ACA", "ACG", "ACT"}},
    {"Phe", {"TTT", "TTC"}},
    {"Ala", {"GCA", "GCC", "GCG", "GCT"}},
    {"Gly", {"GGT", "GGG", "GGA", "GGC"}},
    {"Ile", {"ATC", "ATA", "ATT"}},
    {"Leu", {"TTA", "TTG", "CTC", "CTT", "CTG", "CTA"}},
    {"His", {"CAT", "CAC"}},
    {"Arg", {"CGA", "CGC", "CGG", "CGT", "AGG", "AGA"}},
    {"Trp", {"TGG"}},
    {"Val", {"GTA", "GTC", "GTG", "GTT"}},
    {"Glu", {"GAG", "GAA"}},
    {"Tyr", {"TAT", "TAC"}},
};
