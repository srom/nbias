#include <iostream>
#include <fstream>
#include <cmath>
#include <xtensor/xarray.hpp>
#include "../src/domain/domain.cpp"
#define BOOST_TEST_MODULE "Domain test"
#include <boost/test/unit_test.hpp>

using namespace std;
namespace utf = boost::unit_test;

BOOST_AUTO_TEST_CASE(TestProteinDomains, * utf::tolerance(0.00001))
{
	ifstream f("../fixtures/GCA_000010525.1_pfam.csv");
	ProteinDomains domains(f);

	BOOST_TEST(domains.size() == 2'546);
	BOOST_TEST(domains.Keys().size() == 2'546);

	ProteinDomain domain("PF05175", "MTS");

	auto protein_ids = domains.ProteinIds(domain);

	BOOST_TEST(protein_ids.size() == 11);
	BOOST_TEST(protein_ids[0] == "BAF86903.1");

	string path = "../fixtures/GCA_000010525.1_tri_nucleotide_distance_to_mean.csv";
	ifstream f1(path);
	GeneProbabilies probs(f1, "left");

	xt::xarray<double> probabilities = domains.Probabilities(domain, probs);
	BOOST_TEST(probabilities.size() == 11);
	BOOST_TEST(probabilities[10] == 0.641768);

	xt::xarray<double> rand_probabilities = domains.Probabilities(domain, probs, true);
	BOOST_TEST(rand_probabilities.size() == 11);
	BOOST_TEST(rand_probabilities[10] != probabilities[10]);
}

BOOST_AUTO_TEST_CASE(TestDomainProbability, * utf::tolerance(0.00001))
{
	ProteinDomain domain("PF05175", "MTS", "Methyltransferase small domain");
	DomainProbability domain_probability(domain, log(0.9), log(0.2), 10);

	BOOST_TEST(domain_probability.evidence == 0.653213);
	BOOST_TEST(domain_probability.evidence_strength == "Substantial");

	ProteinDomain domain2("PF00004", "AAA");
	DomainProbability domain_probability2(domain2, log(0.3), log(0.2), 5);
	BOOST_TEST(domain_probability2.evidence_strength == "Weak");

	bool g = domain_probability > domain_probability2;
	BOOST_TEST(g == true);

	auto record_header = DomainProbability::RecordHeader();
	BOOST_TEST(record_header.size() == 8);
	BOOST_TEST(record_header[0] == "id");
	BOOST_TEST(record_header[1] == "query");
	BOOST_TEST(record_header[2] == "description");
	BOOST_TEST(record_header[3] == "log_probability");
	BOOST_TEST(record_header[4] == "log_probability_random");
	BOOST_TEST(record_header[5] == "n_elements");
	BOOST_TEST(record_header[6] == "evidence");
	BOOST_TEST(record_header[7] == "evidence_strength");

	auto record = domain_probability.Record();
	BOOST_TEST(record.size() == 8);
	BOOST_TEST(record[0] == "PF05175");
	BOOST_TEST(record[1] == "MTS");
	BOOST_TEST(record[2] == "Methyltransferase small domain");
	BOOST_TEST(record[3] == "-0.10536052");
	BOOST_TEST(record[4] == "-1.60943791");
	BOOST_TEST(record[5] == "10");
	BOOST_TEST(record[6] == "0.65321251");
	BOOST_TEST(record[7] == "Substantial");

	// Probability <= 0 or > 1: throw runtime error.
	BOOST_CHECK_THROW(DomainProbability(domain, log(0.0), log(0.2), 5), runtime_error);
	BOOST_CHECK_THROW(DomainProbability(domain, log(-0.1), log(0.2), 5), runtime_error);
	BOOST_CHECK_THROW(DomainProbability(domain, log(0.9), log(1.5), 5), runtime_error);
}

BOOST_AUTO_TEST_CASE(TestLoadDomainProbabilities, * utf::tolerance(0.00001)) {
	string path = "../fixtures/GCA_000010525.1_pfam_probability_left.csv";
	vector<DomainProbability> domains = LoadDomainProbabilities(path);

	BOOST_TEST(domains.size() == 2'546);

	size_t ix = 794;
	BOOST_TEST(domains[ix].domain.id == "PF06481");
	BOOST_TEST(domains[ix].domain.query == "COX_ARM");
	BOOST_TEST(domains[ix].domain.description == "COX Aromatic Rich Motif");
	BOOST_TEST(domains[ix].log_probability == log(0.692109));
	BOOST_TEST(domains[ix].log_probability_random == log(0.611623));
	BOOST_TEST(domains[ix].evidence == 0.0536908);
	BOOST_TEST(domains[ix].evidence_strength == "Weak");
}

BOOST_AUTO_TEST_CASE(TestLoadDomainMetadata, * utf::tolerance(0.00001)) {
	string path = "../fixtures/pfam_master.csv";
	auto metadata = LoadDomainMetadata(path);

	BOOST_TEST(metadata.size() == 9);

	auto& [query, description] = metadata["PF00004"];

	BOOST_TEST(query == "AAA");
	BOOST_TEST((
		description == 
		"ATPase family associated with various cellular activities (AAA)"
	));
}
