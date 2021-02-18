#pragma once
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
 * @param rng: A Mersenne Twister pseudo-random generator of 32-bit numbers.
 *
 * @return sequence with synonymous codons randomly swapped.
 */
string SwapSynonymousCodons(const string& sequence, mt19937 rng) {
	if (sequence.size() < 3) {
		return sequence;
	}
	unordered_map<int, uniform_int_distribution<>> uniform_dists;
	string swapped = sequence;
	for (int pos = 0; pos + 3 <= sequence.size(); pos += 3) {
		string codon = swapped.substr(pos, pos + 3);

		if (CODON_TO_AMINO_ACID.find(codon) != CODON_TO_AMINO_ACID.end()) {
			const string& amino_acid = CODON_TO_AMINO_ACID.at(codon);
			const vector<string>& synonymous_codons = AMINO_ACID_TO_CODONS.at(amino_acid);

			auto size = synonymous_codons.size();
			if (size > 1) {
				if (uniform_dists.find(size) == uniform_dists.end()) {
					uniform_dists[size] = uniform_int_distribution<>(0, size - 1);
				}
				auto& distr = uniform_dists[size];
				const string& synonymous_codon = synonymous_codons[distr(rng)];
				swapped.replace(pos, 3, synonymous_codon);
			}
		}
	}
	return swapped;
}

/**
 * Randomly shuffle the order of codons in sequence.
 *
 * @param sequence: coding gene sequence.
 *
 * @param rng: A Mersenne Twister pseudo-random generator of 32-bit numbers.
 *
 * @return sequence with codons randomly shuffled.
 */
string ShuffleCodons(const string& sequence, mt19937 rng) {
	if (sequence.size() < 3) {
		return sequence;
	}
	// Extract codons
	vector<string> sequence_codons;
	for (int pos = 0; pos + 3 <= sequence.size(); pos += 3) {
		sequence_codons.push_back(sequence.substr(pos, 3));
	}

	// Shuffle codons
	shuffle(sequence_codons.begin(), sequence_codons.end(), rng);

	// Reconstruct string
	string shuffled;
	for (const auto& codon : sequence_codons) {
		shuffled += codon;
	}

	// Add characters at the end if necessary
	auto reminder = sequence.size() % 3;
	if (reminder != 0) {
		shuffled += sequence.substr(sequence.size() - reminder);
	}
	return shuffled;
}
