#include <string>
#include <algorithm>

using namespace std;

void ReverseComplement(string& seq) {
	reverse(seq.begin(), seq.end());
	for (size_t i = 0; i < seq.size(); ++i){
		switch (seq[i]) {
			case 'A':
				seq[i] = 'T';
				break;    
			case 'C':
				seq[i] = 'G';
				break;
			case 'G':
				seq[i] = 'C';
				break;
			case 'T':
				seq[i] = 'A';
				break;
			}
	}
}
