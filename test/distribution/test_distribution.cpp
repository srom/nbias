#include <iostream>
#include <fstream>
#include <xtensor/xmath.hpp>
#include "../src/distribution/distribution.cpp"
#define BOOST_TEST_MODULE "Distribution test"
#include <boost/test/unit_test.hpp>

using namespace std;
namespace utf = boost::unit_test;

BOOST_AUTO_TEST_CASE(TestLoadDistributions, * utf::tolerance(0.00001))
{
	ifstream f("../fixtures/distribution.csv");

	Distributions distributions(f);

	BOOST_TEST(distributions.size() == 10);

	auto keys = distributions.Keys();

	BOOST_TEST(keys.size() == 10);
	BOOST_TEST(keys[4] == "GCA_001735525.1");

	auto dist = distributions["GCA_001735525.1"];

	BOOST_TEST(dist.size() == 64);
	BOOST_TEST(dist[1] == 0.018183);
	BOOST_TEST(xt::sum(dist)[0] == 1.0);

	auto values = distributions.Values();
	BOOST_TEST(values(4, 1) == 0.018183);
	BOOST_TEST(xt::sum(xt::row(values, 4))[0] == 1.0);
}
