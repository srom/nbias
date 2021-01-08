#include <iostream>
#include <fstream>
#include <xtensor/xarray.hpp>
#include <xtensor/xrandom.hpp>
#include "../src/domain/domain.cpp"
#define BOOST_TEST_MODULE "Domain test"
#include <boost/test/unit_test.hpp>

using namespace std;
namespace utf = boost::unit_test;

BOOST_AUTO_TEST_CASE(TestProteinDomains, * utf::tolerance(0.00001))
{
	xt::random::seed(444);

	ifstream f("../fixtures/GCA_000010525.1_pfam.csv");
	ProteinDomains domains(f);

	BOOST_TEST(domains.size() == 2'546);
	BOOST_TEST(domains.Keys().size() == 2'546);

	ProteinDomain domain("PF05175.14", "MTS");

	auto protein_ids = domains.ProteinIds(domain);

	BOOST_TEST(protein_ids.size() == 11);
	BOOST_TEST(protein_ids[0] == "BAF86903.1");

	string path = "../fixtures/GCA_000010525.1_tri_nucleotide_distance_to_mean.csv";
	ifstream f1(path);
	GeneProbabilies probs(f1, "left");

	xt::xarray<double> probabilities = domains.Probabilities(domain, probs);
	BOOST_TEST(probabilities.size() == 11);
	BOOST_TEST(probabilities[10] == 0.923240);

	xt::xarray<double> rand_probabilities = domains.Probabilities(domain, probs, true);
	BOOST_TEST(rand_probabilities.size() == 11);
	BOOST_TEST(rand_probabilities[10] != probabilities[10]);
}

BOOST_AUTO_TEST_CASE(TestDomainProbability, * utf::tolerance(0.00001))
{
	ProteinDomain domain("PF05175.14", "MTS");
	DomainProbability domain_probability(domain, 0.9, 0.2);

	BOOST_TEST(domain_probability.evidence == 0.653213);
	BOOST_TEST(domain_probability.evidence_strength == "Substantial");

	ProteinDomain domain2("PF00004.1", "AAA");
	DomainProbability domain_probability2(domain2, 0.3, 0.2);
	BOOST_TEST(domain_probability2.evidence_strength == "Weak");

	bool g = domain_probability > domain_probability2;
	BOOST_TEST(g == true);

	auto record_header = DomainProbability::RecordHeader();
	BOOST_TEST(record_header.size() == 6);
	BOOST_TEST(record_header[0] == "id");
	BOOST_TEST(record_header[1] == "query");
	BOOST_TEST(record_header[2] == "probability");
	BOOST_TEST(record_header[3] == "probability_random");
	BOOST_TEST(record_header[4] == "evidence");
	BOOST_TEST(record_header[5] == "evidence_strength");

	auto record = domain_probability.Record();
	BOOST_TEST(record.size() == 6);
	BOOST_TEST(record[0] == "PF05175.14");
	BOOST_TEST(record[1] == "MTS");
	BOOST_TEST(record[2] == "0.900000");
	BOOST_TEST(record[3] == "0.200000");
	BOOST_TEST(record[4] == "0.653213");
	BOOST_TEST(record[5] == "Substantial");
}
