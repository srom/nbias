#include <iostream>
#include <cmath>
#include <xtensor/xarray.hpp>
#include <xtensor/xmath.hpp>
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

	BOOST_TEST(v == 0.5108256);

	xt::xarray<double> p1 = {
		{0.5, 0.5},
		{(double) 9 / 10, (double) 1 / 10}
	};
	xt::xarray<double> q1  = {
		{(double) 9 / 10, (double) 1 / 10},
		{0.5, 0.5}
	};

	xt::xarray<double> res1 = kl_divergence(p1, q1, 1);

	BOOST_TEST(res1[0] == 0.5108256);
	BOOST_TEST(res1[1] == 0.3680642);

	xt::xarray<double> p0 = {
		{0.5, (double) 9 / 10},
		{0.5, (double) 1 / 10}
	};
	xt::xarray<double> q0  = {
		{(double) 9 / 10, 0.5},
		{(double) 1 / 10, 0.5}
	};
	
	xt::xarray<double> res0 = kl_divergence(p0, q0, 0);

	BOOST_TEST(res0[0] == 0.5108256);
	BOOST_TEST(res0[1] == 0.3680642);

	// Zero value: kl_divergence should throw an exception.
	xt::xarray<double> a = {1.0, 0.0};
	xt::xarray<double> b = {0.5, 0.5};

	BOOST_CHECK_THROW(kl_divergence(a, b), runtime_error);

	xt::xarray<double> w = {
		{0.5, 0.},
		{(double) 9 / 10, (double) 1 / 10}
	};
	xt::xarray<double> x  = {
		{(double) 9 / 10, (double) 1 / 10},
		{0.5, 0.5}
	};

	BOOST_CHECK_THROW(kl_divergence(w, x, 1);, runtime_error);

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

	// Shape mismatch exception
	BOOST_CHECK_THROW(kl_divergence(p, q0);, runtime_error);
	BOOST_CHECK_THROW(kl_divergence(p, q0, 0);, runtime_error);
	BOOST_CHECK_THROW(kl_divergence(p, q0, 1);, runtime_error);
}

BOOST_AUTO_TEST_CASE(TestJensenShannonDistance, * utf::tolerance(0.00001)) {
	xt::xarray<double> p = {0.5, 0.5};
	xt::xarray<double> q = {(double) 9 / 10, (double) 1 / 10};

	double v = jensen_shannon_distance(p, q);

	BOOST_TEST(v == 0.3189815);

	// distance is symmetric
	double v_sym = jensen_shannon_distance(q, p);

	BOOST_TEST(v_sym == 0.3189815);

	xt::xarray<double> p1 = {
		{0.5, 0.5},
		{(double) 9 / 10, (double) 1 / 10}
	};
	xt::xarray<double> q1  = {
		{(double) 9 / 10, (double) 1 / 10},
		{0.5, 0.5}
	};

	xt::xarray<double> res1 = jensen_shannon_distance(p1, q1, 1);

	BOOST_TEST(res1[0] == 0.3189815);
	BOOST_TEST(res1[1] == 0.3189815);

	xt::xarray<double> res1_sym = jensen_shannon_distance(q1, p1, 1);

	BOOST_TEST(res1_sym[0] == 0.3189815);
	BOOST_TEST(res1_sym[1] == 0.3189815);

	// Zero value: jensen_shannon_distance should throw an exception.
	xt::xarray<double> a = {1.0, 0.0};
	xt::xarray<double> b = {0.5, 0.5};

	BOOST_CHECK_THROW(jensen_shannon_distance(a, b), runtime_error);

	// Negative number: jensen_shannon_distance should throw an exception.
	xt::xarray<double> d = {1.0, -0.00001};
	xt::xarray<double> e = {0.5, 0.5};

	BOOST_CHECK_THROW(jensen_shannon_distance(d, e), runtime_error);

	// Shape mismatch exception
	BOOST_CHECK_THROW(jensen_shannon_distance(p, q1);, runtime_error);
	BOOST_CHECK_THROW(jensen_shannon_distance(p, q1, 0);, runtime_error);
	BOOST_CHECK_THROW(jensen_shannon_distance(p, q1, 1);, runtime_error);
}
