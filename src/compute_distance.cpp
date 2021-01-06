/*
Compute Jensen-Shannon distance between an assembly's average tri-nucleotide distribution 
and its genes' tri-nucleotide distributions.
*/
#include <iostream>
#include <string>
#include <boost/program_options.hpp>
#include <xtensor/xarray.hpp>
#include "compute_distance_task.cpp"

using namespace std;
namespace po = boost::program_options;

int main(int ac, char* av[]) {
	try {
		po::options_description desc(
			"Compute Jensen-Shannon distance between an assembly's average tri-nucleotide "
            "distribution and its genes' tri-nucleotide distributions"
		);
		desc.add_options()
            (
            	"help,h", 
            	"Print help message"
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

        const int n_threads = vm["n_threads"].as<int>();
        cerr << "Threads: " << n_threads << endl;

        compute_distance(n_threads);
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
