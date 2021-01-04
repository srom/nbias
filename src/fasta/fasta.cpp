/*
Simple Fasta reader.
*/
#include <fstream>
#include <iostream>
#include <string>
#include <regex>

using namespace std;

const regex _recordStartPattern{R"(^>\s*([^\s]+)\s*(.*)$)"};

struct FastaRecord {
	string id;
	string description;
	string content;

	bool empty() {
		return id.empty() || content.empty();
	}
};

class FastaParser {
	istream& is;
	string id;
	string description;
	string content;
	
	public:
		FastaParser(istream& s) :is{s} {}

		bool Get(FastaRecord& record) {
			string line;
			while (getline(is, line)) {
				if (line.empty()) {
					continue;
				}
				if (line[0] == '>') {
					if (!id.empty() && !content.empty()) {
						record = FastaRecord{id, description, content};
						id.clear();
						description.clear();
						content.clear();
					}
					smatch m;
					if (!regex_match(line, m, _recordStartPattern)) {
						throw runtime_error("Invalid record metadata: " + line);
					}
					id = m[1];
					if (m.size() > 2) {
						description = m[2];
					}
					if (!record.empty()) {
						return true;
					}
				} else if (!id.empty()) {
					content += line;
				}
			}
			if (!id.empty() && !content.empty()) {
				record = FastaRecord{id, description, content};
				id.clear();
				description.clear();
				content.clear();
				return true;
			} else {
				record = FastaRecord{};
			}
			return false;
		}

		vector<FastaRecord> GetAll() {
			vector<FastaRecord> records;
			FastaRecord record{};
			while (Get(record)) {
				records.push_back(record);
			};
			return records;
		}
};
