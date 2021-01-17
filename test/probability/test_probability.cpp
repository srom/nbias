#include <iostream>
#include <fstream>
#include <cmath>
#include <xtensor/xarray.hpp>
#include <xtensor/xmath.hpp>
#include "../src/probability/probability.cpp"
#include "../src/probability/probability_util.cpp"
#define BOOST_TEST_MODULE "Probability test"
#include <boost/test/unit_test.hpp>

using namespace std;
namespace utf = boost::unit_test;

BOOST_AUTO_TEST_CASE(TestProductRule, * utf::tolerance(0.00001))
{
	xt::xarray<double> p = {0.9, 0.9, 0.1};
	BOOST_TEST(product_rule(p) == 0.081);

	xt::xarray<double> p_matrix = {
		{0.9, 0.9, 0.1},
		{0.4, 0.2, 0.1},
	};

	xt::xarray<double> r0 = product_rule(p_matrix, 0);
	xt::xarray<double> r1 = product_rule(p_matrix, 1);

	BOOST_TEST(r0.size() == 3);
	BOOST_TEST(r1.size() == 2);

	BOOST_TEST(r0[0] == 0.36);
	BOOST_TEST(r0[1] == 0.18);
	BOOST_TEST(r0[2] == 0.01);

	BOOST_TEST(r1[0] == 0.081);
	BOOST_TEST(r1[1] == 0.008);

	// One element is zero: product is zero.
	xt::xarray<double> q = {0., 0.9, 0.1};
	BOOST_TEST(product_rule(q) == 0);

	// One element is negative: product is nan.
	xt::xarray<double> r = {-0.001, 0.9, 0.1};
	BOOST_TEST(isnan(product_rule(r)));
}

BOOST_AUTO_TEST_CASE(TestProductRuleLog, * utf::tolerance(0.00001))
{
	xt::xarray<double> p = {0.9, 0.9, 0.1};
	xt::xarray<double> log_p = xt::log(p);

	BOOST_TEST(product_rule_log(log_p) == log(0.081));

	xt::xarray<double> p_matrix = {
		{0.9, 0.9, 0.1},
		{0.4, 0.2, 0.1},
	};
	xt::xarray<double> log_p_matrix = xt::log(p_matrix);

	xt::xarray<double> r1 = product_rule_log(log_p_matrix, 1);

	BOOST_TEST(r1.size() == 2);
	BOOST_TEST(r1[0] == log(0.081));
	BOOST_TEST(r1[1] == log(0.008));

	// One element is zero: log product is -inf.
	xt::xarray<double> q = {0., 0.9, 0.1};
	xt::xarray<double> log_q = xt::log(q);
	BOOST_TEST(isinf(product_rule_log(log_q)));

	// One element is negative: product is nan.
	xt::xarray<double> r = {-0.001, 0.9, 0.1};
	xt::xarray<double> log_r = xt::log(r);
	BOOST_TEST(isnan(product_rule_log(log_r)));
}

BOOST_AUTO_TEST_CASE(TestProductRuleNormalized, * utf::tolerance(0.00001))
{
	xt::xarray<double> p = {0.9, 0.9, 0.1};
	BOOST_TEST(product_rule_normalized(p) == 0.432674);

	xt::xarray<double> p_matrix = {
		{0.9, 0.9, 0.1},
		{0.4, 0.2, 0.1},
	};

	xt::xarray<double> r0 = product_rule_normalized(p_matrix, 0);
	xt::xarray<double> r1 = product_rule_normalized(p_matrix, 1);

	BOOST_TEST(r0.size() == 3);
	BOOST_TEST(r1.size() == 2);

	BOOST_TEST(r0[0] == 0.600000);
	BOOST_TEST(r0[1] == 0.424264);
	BOOST_TEST(r0[2] == 0.100000);

	BOOST_TEST(r1[0] == 0.432674);
	BOOST_TEST(r1[1] == 0.200000);

	// One element is zero: product is zero.
	xt::xarray<double> q = {0., 0.9, 0.1};
	BOOST_TEST(product_rule_normalized(q) == 0);

	// One element is negative: product is nan.
	xt::xarray<double> r = {-0.001, 0.9, 0.1};
	BOOST_TEST(isnan(product_rule_normalized(r)));
}

BOOST_AUTO_TEST_CASE(TestMarginalization, * utf::tolerance(0.00001))
{
	xt::xarray<double> prior = {0.1, 0.5, 0.4};
	xt::xarray<double> likelihood = {0.9, 0.9, 0.1};
	BOOST_TEST(marginalization(prior, likelihood) == 0.580000);

	xt::xarray<double> prior_mat = {
		{0.100, 0.500, 0.400},
		{0.333, 0.333, 0.333},
	};
	xt::xarray<double> likelihood_mat = {
		{0.9, 0.9, 0.1},
		{0.4, 0.2, 0.1},
	};

	xt::xarray<double> r = marginalization(prior_mat, likelihood_mat, 1);

	BOOST_TEST(r.size() == 2);
	BOOST_TEST(r[0] == 0.580000);
	BOOST_TEST(r[1] == 0.233100);

	// One element is zero: should be handled transparently.
	BOOST_TEST(marginalization({0.5, 0.5}, {0.0, 0.9}) == 0.45);

	// One element is negative: marginalization is nan.
	BOOST_TEST(isnan(marginalization({0.5, 0.5}, {-0.001, 0.9})));
}

BOOST_AUTO_TEST_CASE(TestMarginalizationLog, * utf::tolerance(0.00001))
{
	xt::xarray<double> prior = {0.1, 0.5, 0.4};
	xt::xarray<double> likelihood = {0.9, 0.9, 0.1};
	xt::xarray<double> log_prior = xt::log(prior);
	xt::xarray<double> log_likelihood = xt::log(likelihood);

	BOOST_TEST(marginalization_log(log_prior, log_likelihood) == log(0.58));

	xt::xarray<double> prior_mat = {
		{0.100, 0.500, 0.400},
		{0.333, 0.333, 0.333},
	};
	xt::xarray<double> likelihood_mat = {
		{0.9, 0.9, 0.1},
		{0.4, 0.2, 0.1},
	};
	xt::xarray<double> log_prior_mat = xt::log(prior_mat);
	xt::xarray<double> log_likelihood_mat = xt::log(likelihood_mat);

	xt::xarray<double> r = marginalization_log(log_prior_mat, log_likelihood_mat, 1);

	BOOST_TEST(r.size() == 2);
	BOOST_TEST(r[0] == log(0.5800));
	BOOST_TEST(r[1] == log(0.2331));

	xt::xarray<double> p = {0.5, 0.5};
	xt::xarray<double> l = {0.0, 0.9};
	xt::xarray<double> log_p = xt::log(p);
	xt::xarray<double> log_l = xt::log(l);

	// One element is zero: should be handled transparently.
	BOOST_TEST(marginalization_log(log_p, log_l) == log(0.45));

	xt::xarray<double> l_neg = {-0.001, 0.9};
	xt::xarray<double> log_l_neg = xt::log(l_neg);

	// One element is negative: marginalization is nan.
	BOOST_TEST(isnan(marginalization_log(log_p, log_l_neg)));
}

BOOST_AUTO_TEST_CASE(TestComputeProbabilityFromDistance, * utf::tolerance(0.00001))
{
	vector<double> distances{0.95, 0.92, 0.71};
	xt::xarray<double> p = compute_probability_from_distance(distances, true);
	xt::xarray<double> q = compute_probability_from_distance(distances, false);

	BOOST_TEST(p.size() == 3);
	BOOST_TEST(q.size() == 3);

	BOOST_TEST(p[0] == 0.05);
	BOOST_TEST(p[1] == 0.08);
	BOOST_TEST(p[2] == 0.29);

	BOOST_TEST(p[0] == 1.0 - q[0]);
	BOOST_TEST(p[1] == 1.0 - q[1]);
	BOOST_TEST(p[2] == 1.0 - q[2]);
}

BOOST_AUTO_TEST_CASE(TestGeneProbabilities, * utf::tolerance(0.00001))
{
	string path = "../fixtures/GCA_000010525.1_tri_nucleotide_distance_to_mean.csv";
	ifstream f1(path);
	ifstream f2(path);

	GeneProbabilies p(f1, "left");
	GeneProbabilies q(f2, "right");

	BOOST_TEST(p["BAF86014.1"] == 1.0 - 0.080520);
	BOOST_TEST(q["BAF86014.1"] == 0.080520);

	BOOST_TEST(q["BAF86014.1"] == q.Get("BAF86014.1"));
	BOOST_TEST(q["BAF86014.1"] != q.Get("BAF86014.1", true));

	BOOST_TEST(p.size() == 4'717);
	BOOST_TEST(p.Values().size() == 4'717);

	BOOST_TEST(p.size() == 4'717);
	BOOST_TEST(q.Values().size() == 4'717);
	BOOST_TEST(q.RandomValues().size() == 4'717);

	BOOST_TEST(p.Values()[4'716] == 1.0 - 0.074669);
	BOOST_TEST(q.Values()[4'716] == 0.074669);

	BOOST_TEST(p.Keys().size() == 4'717);
	BOOST_TEST(q.Keys().size() == 4'717);

	BOOST_TEST(p.IndexOf("BAF86000.1") == 1);
}

BOOST_AUTO_TEST_CASE(TestMakeUniformPrior, * utf::tolerance(0.00001))
{
	xt::xarray<double> p = make_uniform_prior(4);

	BOOST_TEST(p.size() == 4);
	BOOST_TEST(p[0] == 0.25);
	BOOST_TEST(p[1] == 0.25);
	BOOST_TEST(p[2] == 0.25);
	BOOST_TEST(p[3] == 0.25);

	xt::xarray<double> q = make_uniform_prior(6, 4);
	
	BOOST_TEST(q.size() == 6 * 4);
	BOOST_TEST(q.shape()[0] == 6);
	BOOST_TEST(q.shape()[1] == 4);
	BOOST_TEST(q(0, 2) == 0.25);
	BOOST_TEST(q(1, 3) == 0.25);
}
