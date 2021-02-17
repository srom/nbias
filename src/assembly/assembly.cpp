/*
Assembly objects represent NCBI genome assemblies.
Assemblies hold Assembly objects and provide methods for easy access.
*/
#include <iostream>
#include <string>
#include <unordered_map>
#include "csv.hpp"

using namespace std;
using namespace csv;

struct Assembly {
	string assembly_accession;
	int taxid;
	int species_taxid;
	string organism_name;
	string domain;
	string phylum;
	string classification;
	string order;
	string family;
	string genus;
	string species;
	string strain;
	string assembly_level;
};

class Assemblies {
	vector<string> assembly_accessions;
	unordered_map<string, Assembly> data;

	public:
		Assemblies(string path, bool complete_genome_only = false) :assembly_accessions{}, data{} {
			CSVReader reader(path);
			for (CSVRow& row: reader) {
				string phylum = row["phylum"].get<string>();
				string assembly_level = row["assembly_level"].get<string>();
				string assembly_accession = row["assembly_accession"].get<string>();

				bool not_a_complete_genome = (
					assembly_level != "Complete Genome" or
					phylum.size() == 0
				);
				if (complete_genome_only and not_a_complete_genome) {
					continue;
				}

				Assembly assembly{
					assembly_accession,
					row["taxid"].get<int>(),
					row["species_taxid"].get<int>(),
					row["organism_name"].get<string>(),
					row["domain"].get<string>(),
					phylum,
					row["class"].get<string>(),
					row["order"].get<string>(),
					row["family"].get<string>(),
					row["genus"].get<string>(),
					row["species"].get<string>(),
					row["strain"].get<string>(),
					assembly_level,
				};
				assembly_accessions.push_back(assembly_accession);
				data[assembly_accession] = assembly;
			}
		}

		vector<string> GetIds() {
			return assembly_accessions;
		}

		bool Has(const string assembly_accession) {
			return data.find(assembly_accession) != data.end();
		}

		Assembly& Get(const string assembly_accession) {
			return data[assembly_accession];
		}

		size_t Size() {
			return assembly_accessions.size();
		}
};
