#include <string>
#include "../src/task/context.cpp"
#define BOOST_TEST_MODULE "Context test"
#include <boost/test/unit_test.hpp>

using namespace std;
namespace utf = boost::unit_test;

BOOST_AUTO_TEST_CASE(TestLoadDistributions, * utf::tolerance(0.00001))
{
	const string path = "../fixtures/domain_prob_conf.json";

	DomainProbabilityContext ctx = parse_domain_probability_context(path);

	BOOST_TEST(ctx.kind == "amino-acid");
	BOOST_TEST(ctx.query == "pfam");
	BOOST_TEST(ctx.tail == "left");
	BOOST_TEST(ctx.n_threads == 5);
	BOOST_TEST(ctx.complete_genome_only == true);
	BOOST_TEST(ctx.distance_to_mean_suffix == "_amino_acid_distance_to_mean.csv");
	BOOST_TEST(ctx.assembly_output_folder == "confounders/length/assembly");
	BOOST_TEST(ctx.phylum_output_folder == "confounders/length/phylum");
	BOOST_TEST(ctx.superkingdom_output_folder == "confounders/length/superkingdom");
	BOOST_TEST(ctx.overall_output_folder == "confounders/length");

	DomainProbabilityContext ctx2 = get_context_from_parameters(
		"amino-acid",
		"tigr",
		"right",
		3
	);

	BOOST_TEST(ctx2.kind == "amino-acid");
	BOOST_TEST(ctx2.query == "tigr");
	BOOST_TEST(ctx2.tail == "right");
	BOOST_TEST(ctx2.n_threads == 3);
	BOOST_TEST(ctx2.complete_genome_only == false);
	BOOST_TEST(ctx2.distance_to_mean_suffix == "_amino_acid_distance_to_mean.csv");
	BOOST_TEST(ctx2.assembly_output_folder == "");
	BOOST_TEST(ctx2.phylum_output_folder == "");
	BOOST_TEST(ctx2.superkingdom_output_folder == "");
	BOOST_TEST(ctx2.overall_output_folder == "");
}
