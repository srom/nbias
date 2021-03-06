#include <iostream>
#include <cmath>
#include <xtensor/xarray.hpp>
#include <xtensor/xmath.hpp>
#include <xtensor/xview.hpp>
#include "../src/distance/distance.cpp"
#define BOOST_TEST_MODULE "Distance test"
#include <boost/test/unit_test.hpp>

using namespace std;
namespace utf = boost::unit_test;

BOOST_AUTO_TEST_CASE(TestKLDivergence, * utf::tolerance(0.00001))
{
	xt::xarray<double> p = {0.5, 0.5};
	xt::xarray<double> q = {(double) 9 / 10, (double) 1 / 10};

	double v = kl_divergence(p, q);

	BOOST_TEST(v == 0.7369656);

	double v2 = kl_divergence(q, p);

	BOOST_TEST(v2 == 0.5310044);

	xt::xarray<double> p1 = {
		{0.5, 0.5},
		{(double) 9 / 10, (double) 1 / 10}
	};
	xt::xarray<double> q1  = {
		{(double) 9 / 10, (double) 1 / 10},
		{0.5, 0.5}
	};

	xt::xarray<double> res1 = kl_divergence(p1, q1, 1);

	BOOST_TEST(res1[0] == 0.7369656);
	BOOST_TEST(res1[1] == 0.5310044);

	xt::xarray<double> p0 = {
		{0.5, (double) 9 / 10},
		{0.5, (double) 1 / 10}
	};
	xt::xarray<double> q0  = {
		{(double) 9 / 10, 0.5},
		{(double) 1 / 10, 0.5}
	};
	
	xt::xarray<double> res0 = kl_divergence(p0, q0, 0);

	BOOST_TEST(res0[0] == 0.7369656);
	BOOST_TEST(res0[1] == 0.5310044);

	xt::xarray<double> a = {1.0, 0.0};
	xt::xarray<double> b = {0.5, 0.5};

	BOOST_TEST(kl_divergence(a, b) == 1.0);

	// Throw: values in q are not allowed to be zero
	BOOST_CHECK_THROW(kl_divergence(b, a), runtime_error);

	// Negative number: kl_divergence should throw an exception.
	xt::xarray<double> d = {1.0, -0.00001};
	xt::xarray<double> e = {0.5, 0.5};

	BOOST_CHECK_THROW(kl_divergence(d, e), runtime_error);

	xt::xarray<double> y = {
		{-0.001, (double) 9 / 10},
		{0.5, (double) 1 / 10}
	};
	xt::xarray<double> z  = {
		{(double) 9 / 10, 0.5},
		{(double) 1 / 10, 0.5}
	};
	
	BOOST_CHECK_THROW(kl_divergence(y, z, 0);, runtime_error);

	// All values are zero: throw
	xt::xarray<double> pzero = {0.0, 0.0};

	BOOST_CHECK_THROW(kl_divergence(pzero, q), runtime_error);
	BOOST_CHECK_THROW(kl_divergence(q, pzero), runtime_error);

	xt::xarray<double> p1zero = {{0.0, 0.0}, {0.5, 0.5}};

	BOOST_CHECK_THROW(kl_divergence(p1zero, q1, 1), runtime_error);

	// Values do not sum to 1: throw
	xt::xarray<double> p_not_one = {0.33, 0.1};
	
	BOOST_CHECK_THROW(kl_divergence(p_not_one, q), runtime_error);
	BOOST_CHECK_THROW(kl_divergence(q, p_not_one), runtime_error);

	xt::xarray<double> p_not_one_1 = {{0.33, 0.1}, {0.5, 0.5}};
	
	BOOST_CHECK_THROW(kl_divergence(p_not_one_1, q1, 1), runtime_error);

	// Shape mismatch exception
	BOOST_CHECK_THROW(kl_divergence(p, q0);, runtime_error);
	BOOST_CHECK_THROW(kl_divergence(p, q0, 0);, runtime_error);
	BOOST_CHECK_THROW(kl_divergence(p, q0, 1);, runtime_error);
}

BOOST_AUTO_TEST_CASE(TestJensenShannonDistance, * utf::tolerance(0.00001)) {
	xt::xarray<double> p = {0.5, 0.5};
	xt::xarray<double> q = {(double) 9 / 10, (double) 1 / 10};

	double v = jensen_shannon_distance(p, q);

	BOOST_TEST(v == 0.3831359);

	// distance is symmetric
	double v_sym = jensen_shannon_distance(q, p);

	BOOST_TEST(v_sym == 0.3831359);

	xt::xarray<double> p1 = {
		{0.5, 0.5},
		{(double) 9 / 10, (double) 1 / 10}
	};
	xt::xarray<double> q1  = {
		{(double) 9 / 10, (double) 1 / 10},
		{0.5, 0.5}
	};

	xt::xarray<double> res1 = jensen_shannon_distance(p1, q1, 1);

	BOOST_TEST(res1[0] == 0.3831359);
	BOOST_TEST(res1[1] == 0.3831359);

	xt::xarray<double> res1_sym = jensen_shannon_distance(q1, p1, 1);

	BOOST_TEST(res1_sym[0] == 0.3831359);
	BOOST_TEST(res1_sym[1] == 0.3831359);

	xt::xarray<double> a = {1.0, 0.0};
	xt::xarray<double> b = {0.5, 0.5};

	BOOST_TEST(jensen_shannon_distance(a, b) == 0.5579230);
	BOOST_TEST(jensen_shannon_distance(b, a) == 0.5579230);

	xt::xarray<double> aa = {1.0, 0.0};
	xt::xarray<double> bb = {0.0, 1.0};

	BOOST_TEST(jensen_shannon_distance(aa, bb) == 1.0);

	// All values are zero: throw
	xt::xarray<double> pzero = {0.0, 0.0};

	BOOST_CHECK_THROW(jensen_shannon_distance(pzero, q), runtime_error);
	BOOST_CHECK_THROW(jensen_shannon_distance(q, pzero), runtime_error);

	xt::xarray<double> p1zero = {{0.0, 0.0}, {0.5, 0.5}};

	BOOST_CHECK_THROW(jensen_shannon_distance(p1zero, q1, 1), runtime_error);

	// Values do not sum to 1: throw
	xt::xarray<double> p_not_one = {0.33, 0.1};
	
	BOOST_CHECK_THROW(jensen_shannon_distance(p_not_one, q), runtime_error);
	BOOST_CHECK_THROW(jensen_shannon_distance(q, p_not_one), runtime_error);

	xt::xarray<double> p_not_one_1 = {{0.33, 0.1}, {0.5, 0.5}};
	
	BOOST_CHECK_THROW(jensen_shannon_distance(p_not_one_1, q1, 1), runtime_error);

	// Negative number: jensen_shannon_distance should throw an exception.
	xt::xarray<double> d = {1.0, -0.00001};
	xt::xarray<double> e = {0.5, 0.5};

	BOOST_CHECK_THROW(jensen_shannon_distance(d, e), runtime_error);

	// Shape mismatch exception
	BOOST_CHECK_THROW(jensen_shannon_distance(p, q1);, runtime_error);
	BOOST_CHECK_THROW(jensen_shannon_distance(p, q1, 0);, runtime_error);
	BOOST_CHECK_THROW(jensen_shannon_distance(p, q1, 1);, runtime_error);
}

BOOST_AUTO_TEST_CASE(TestDuplicateRow, * utf::tolerance(0.00001)) {
	xt::xarray<double> a = {0.5, 0.5};

	auto b = duplicate_rows(a, 10);

	BOOST_TEST(b.shape()[0] == 10);
	BOOST_TEST(xt::row(b, 0)[0] == a[0]);
	BOOST_TEST(xt::row(b, 0)[1] == a[1]);
	BOOST_TEST(xt::row(b, 5)[0] == a[0]);
	BOOST_TEST(xt::row(b, 5)[1] == a[1]);
}
