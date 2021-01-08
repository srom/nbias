/*
Compute probability & bayesian evidence that a Pfam or TIGR protein domain is near 
their respective average tri-nucleotide distribution.

Conversely, the same quantities an be computed for domains consistently far from the mean.
*/
#include <iostream>
#include <string>
#include <algorithm>
#include <boost/program_options.hpp>
#include "domain_probability_task.cpp"

using namespace std;
namespace po = boost::program_options;

int main(int ac, char* av[]) {
	try {
		po::options_description desc(
			"Compute probability & bayesian evidence that a Pfam or TIGR "
            "protein domain is near (or far) from their respective "
            "average tri-nucleotide distribution"
		);
		desc.add_options()
            (
            	"help,h", 
            	"Print help message"
            )
            (
                "query,q", 
                po::value<string>()->default_value(""), 
                "Protein domain source = Pfam or TIGR (case insensitive)"
            )
            (
                "tail,d", 
                po::value<string>()->default_value("left"), 
                "Tail of the distribution of distances to tri-nucleotide mean on which to focus: "
                "left (closest to min distance 0) or right (closest to max distance 1)"
            )
            (
            	"n_threads,t", 
            	po::value<int>()->default_value(4),
            	"Number of threads to use."
            )
            (
                "seed,s", 
                po::value<int>()->default_value(42),
                "Random seed to use."
            )
        ;

        po::variables_map vm;
        po::store(po::parse_command_line(ac, av, desc), vm);
        po::notify(vm);

        if (vm.count("help")) {
            cerr << desc << endl;
            return 0;
        }

        auto query = vm["query"].as<string>();
        transform(query.begin(), query.end(), query.begin(), ::tolower);

        cerr << "Query: " << query << endl;
        if (query != "pfam" && query != "tigr") {
            cerr << "Error: --query must be one of pfam, tigr" << endl;
            cerr << desc << endl;
            return 1;
        }

        auto tail = vm["tail"].as<string>();
        transform(tail.begin(), tail.end(), tail.begin(), ::tolower);
        cerr << "Tail: " << tail << endl;
        if (tail != "left" && tail != "right") {
            cerr << "Error: --tail must be one of left, right" << endl;
            cerr << desc << endl;
            return 1;
        }

        const int n_threads = vm["n_threads"].as<int>();
        cerr << "Threads: " << n_threads << endl;

        const int seed = vm["seed"].as<int>();
        cerr << "Random seed: " << seed << endl;

        compute_domain_probabilities(query, tail, n_threads, seed);
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
