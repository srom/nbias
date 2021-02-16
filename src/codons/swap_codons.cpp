#include <numeric>
#include <string>
#include <random>
#include <algorithm>
#include <iostream>
#include "./codons.cpp"

using namespace std;

/**
 * Randomly swap synonymous codons in sequence.
 *
 * @param sequence: coding gene sequence.
 *
 * @param seed: Random seed. Randomly set if negative or not set.
 *
 * @return sequence with synonymous codons randomly swapped.
 */
string SwapSynonymousCodons(const string& sequence, const int seed = -1) {
	if (sequence.size() < 3) {
		return sequence;
	}
	// Set random seed
	random_device rd;
	mt19937 gen(rd());
	if (seed >= 0) {
		gen.seed(seed);
	}
	string swapped = sequence;
	for (int pos = 0; pos + 3 <= sequence.size(); pos += 3) {
		string codon = swapped.substr(pos, pos + 3);

		if (CODON_TO_AMINO_ACID.find(codon) != CODON_TO_AMINO_ACID.end()) {
			const string& amino_acid = CODON_TO_AMINO_ACID.at(codon);
			const vector<string>& synonymous_codons = AMINO_ACID_TO_CODONS.at(amino_acid);

			if (synonymous_codons.size() > 1) {
				uniform_int_distribution<> distr(0, synonymous_codons.size() - 1);
				const string& synonymous_codon = synonymous_codons[distr(gen)];
				swapped.replace(pos, 3, synonymous_codon);
			}
		}
		pos += 3;
	}
	return swapped;
}

/**
 * Randomly shuffle the order of codons in sequence.
 *
 * @param sequence: coding gene sequence.
 *
 * @param seed: Random seed. Randomly set if negative or not set.
 *
 * @return sequence with codons randomly shuffled.
 */
string ShuffleCodons(const string& sequence, const int seed = -1) {
	if (sequence.size() < 3) {
		return sequence;
	}
	// Set random seed
	random_device rd;
	mt19937 gen(rd());
	if (seed >= 0) {
		gen.seed(seed);
	}

	// Extract codons
	vector<string> sequence_codons;
	for (int pos = 0; pos + 3 <= sequence.size(); pos += 3) {
		sequence_codons.push_back(sequence.substr(pos, 3));
	}

	// Shuffle codons
	shuffle(sequence_codons.begin(), sequence_codons.end(), gen);

	// Reconstruct string
	string shuffled = accumulate(sequence_codons.begin(), sequence_codons.end(), ""s);

	// Add characters at the end if necessary
	if (sequence.size() % 3 != 0) {
		auto reminder = sequence.size() % 3;
		shuffled += sequence.substr(sequence.size() - reminder);
	}
	return shuffled;
}
