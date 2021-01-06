/*
Compute tri-nucleotide bias for assemblies in data/sequences. 
*/
#include <iostream>
#include <string>
#include <boost/program_options.hpp>
#include "tri_bias_task.cpp"

using namespace std;
namespace po = boost::program_options;

int main(int ac, char* av[]) {
	try {
		po::options_description desc(
			"Compute tri-nucleotide distribution of sequences in data/sequences"
		);
		desc.add_options()
            (
            	"help,h", 
            	"Print help message"
            )
            (
            	"level,l", 
            	po::value<string>()->default_value(""), 
            	"One of \"genome\" (compute genome-wide), "
            	"\"genes\" (compute on genes only), "
            	"\"cds\" (compute distribution for each individual CDS)."
            )
            (
            	"overlap,o", 
            	po::bool_switch(), 
            	"Count kmers with overlap"
            	"(i.e. consider all reading frames)"
            )
            (
            	"n_threads,t", 
            	po::value<int>()->default_value(4),
            	"Number of threads to use."
            )
        ;

        po::variables_map vm;
        po::store(po::parse_command_line(ac, av, desc), vm);
        po::notify(vm);

        if (vm.count("help")) {
            cerr << desc << endl;
            return 0;
        }

        auto level = vm["level"].as<string>();
        if (level == "genome" || level == "genes" || level == "cds") {
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

        const int n_threads = vm["n_threads"].as<int>();
        cerr << "Threads: " << n_threads << endl;

        if (level == "genome") {
        	run_assembly_count(true, overlap, n_threads);
        } else if (level == "genes") {
        	run_assembly_count(false, overlap, n_threads);
        } else if (level == "cds") {
        	run_cds_count(overlap, n_threads);
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
