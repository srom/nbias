/*
Compute tri-nucleotide or amino acid distributions 
of assemblies from the data/sequences folder.
*/
#include <iostream>
#include <string>
#include <boost/program_options.hpp>
#include "./task/compute_distribution_task.cpp"

using namespace std;
namespace po = boost::program_options;

int main(int ac, char* av[]) {
	try {
		po::options_description desc(
			"Compute tri-nucleotide or amino acid distribution "
            "of sequences in folder data/sequences"
		);
		desc.add_options()
            (
            	"help,h", 
            	"Print help message"
            )
            (
                "kind,k", 
                po::value<string>()->default_value(""), 
                "One of \"tri-nucleotide\" or \"amino-acid\""
            )
            (
            	"level,l", 
            	po::value<string>()->default_value(""), 
            	"One of \"genome\" (compute distribution genome-wide), "
            	"\"cds\" (compute distribution for each individual coding sequences)"
            )
            (
            	"overlap,o", 
            	po::bool_switch(), 
            	"Count kmers with overlap "
            	"(i.e. consider all reading frames)"
            )
            (
                "reverse_complement,r", 
                po::bool_switch(), 
                "Consider reverse complement as well (only relevant for kind = tri-nucleotide)"
            )
            (
            	"n_threads,t", 
            	po::value<int>()->default_value(4),
            	"Number of threads to use"
            )
        ;

        po::variables_map vm;
        po::store(po::parse_command_line(ac, av, desc), vm);
        po::notify(vm);

        if (vm.count("help")) {
            cerr << desc << endl;
            return 0;
        }

        auto kind = vm["kind"].as<string>();
        if (kind == "tri-nucleotide" || kind == "amino-acid") {
            cerr << "Kind: " << kind << endl;
        } else if (kind.empty()) {
            cerr << "Error: parameter \"kind\" not set.";
            cerr << " See --help for usage." << endl;
            return 1;
        } else {
            cerr << "Error: unknown value for parameter \"kind\": \"" << kind;
            cerr << "\". See --help for usage." << endl;
            return 1;
        }

        auto level = vm["level"].as<string>();
        if (level == "genome" || level == "cds") {
        	cerr << "Level: " << level << endl;
        } else if (level.empty()) {
        	cerr << "Error: parameter \"level\" not set.";
            cerr << " See --help for usage." << endl;
            return 1;
        } else {
        	cerr << "Error: unknown value for parameter \"level\": \"" << level;
        	cerr << "\". See --help for usage." << endl;
        	return 1;
        }

        const bool overlap = vm["overlap"].as<bool>();
        if (overlap) {
            cerr << "Overlap: true" << endl;
        } else {
            cerr << "Overlap: false" << endl;
        }

        const bool reverse_complement = vm["reverse_complement"].as<bool>();
        if (reverse_complement) {
            cerr << "Reverse complement: true" << endl;
        } else {
            cerr << "Reverse complement: false" << endl;
        }

        if (reverse_complement && kind == "amino-acid") {
            cerr << "Error: Cannot use reverse_complement option with kind = amino-acid";
            return 1;
        }

        const int n_threads = vm["n_threads"].as<int>();
        cerr << "Threads: " << n_threads << endl;

        if (level == "genome") {
        	run_genome_wide_count(kind, overlap, reverse_complement, n_threads);
        } else if (level == "cds") {
        	run_cds_count(kind, overlap, reverse_complement, n_threads);
        }
	}
	catch (exception& e) {
		cerr << "Exception raised: " << e.what() << endl;
		return 1;
	}
	catch (...) {
		cerr << "Unknown exception raised" << endl;
		return 1;
	}
	return 0;
}
