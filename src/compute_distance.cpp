/*
Compute Jensen-Shannon distance between an assembly's average tri-nucleotide 
or amino acid distribution and its genes' tri-nucleotide or amino acid distributions.
*/
#include <iostream>
#include <string>
#include <boost/program_options.hpp>
#include "./task/compute_distance_task.cpp"

using namespace std;
namespace po = boost::program_options;

int main(int ac, char* av[]) {
	try {
		po::options_description desc(
			"Compute Jensen-Shannon distance between an assembly's average "
            "tri-nucleotide (or amino acid) distribution and its genes' "
            "tri-nucleotide (or amino acid) distributions"
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

        cerr << "Computing Jensen-Shannon distance" << endl;

        string kind = vm["kind"].as<string>();
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

        const int n_threads = vm["n_threads"].as<int>();
        cerr << "Threads: " << n_threads << endl;

        compute_distance(kind, n_threads);
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
