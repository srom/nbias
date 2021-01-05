/*
Implementation of KL divergence and Jensen-Shannon distance.
*/
#include <iostream>
#include <algorithm>
#include <cmath>
#include <xtensor/xarray.hpp>
#include <xtensor/xmath.hpp>
#include <xtensor/xcontainer.hpp>

using namespace std;

void kl_divergence_argument_check(
	const xt::xarray<double>& p, 
	const xt::xarray<double>& q
) {
	if (p.shape() != q.shape()) {
		throw runtime_error("KL divergence: shape mismatch");
	}
	auto checkNeg = [](auto&& a) { return a < (double) 0; };
	if (find_if(p.begin(), p.end(), checkNeg) != p.end()) {
		throw runtime_error("KL divergence: element in p is < 0");
	}
	auto checkNegEq = [](auto&& a) { return a <= (double) 0; };
	if (find_if(q.begin(), q.end(), checkNegEq) != q.end()) {
		throw runtime_error("KL divergence: element in q is <= 0");
	}
	double tol = 0.001;
	size_t expected = p.shape()[0];
	if (p.shape().size() == 1) {
		expected = 1;
	}
	auto s1 = xt::sum(p)[0];
	if (fabs(s1 - expected) > tol) {
		throw runtime_error("KL divergence: values in p do not sum to 1");
	}
	auto s2 = xt::sum(q)[0];
	if (fabs(s2 - expected) > tol) {
		throw runtime_error("KL divergence: values in q do not sum to 1");
	}
}

double _kl_divergence_double_ec(
	const xt::xarray<double>& p, 
	const xt::xarray<double>& q,
	const bool error_check
) {
	if (error_check) {
		kl_divergence_argument_check(p, q);
	}
	return xt::sum(p * (xt::nan_to_num(xt::log(p)) - xt::log(q)))[0];
}

double kl_divergence(
	const xt::xarray<double>& p, 
	const xt::xarray<double>& q
) {
	return _kl_divergence_double_ec(p, q, true);
}

xt::xarray<double> _kl_divergence_array_ec(
	const xt::xarray<double>& p, 
	const xt::xarray<double>& q,
	const size_t axis,
	const bool error_check
) {
	if (error_check) {
		kl_divergence_argument_check(p, q);
	}
	return xt::eval(xt::sum(p * (xt::nan_to_num(xt::log(p)) - xt::log(q)), move(axis)));
}

xt::xarray<double> kl_divergence(
	const xt::xarray<double>& p, 
	const xt::xarray<double>& q,
	const size_t axis
) {
	return _kl_divergence_array_ec(p, q, axis, true);
}

void jensen_shannon_distance_argument_check(
	const xt::xarray<double>& p, 
	const xt::xarray<double>& q
) {
	if (p.shape() != q.shape()) {
		throw runtime_error("Jensen-Shannon distance: shape mismatch");
	}
	auto checkNeg = [](auto&& a) { return a < (double) 0; };
	if (find_if(p.begin(), p.end(), checkNeg) != p.end()) {
		throw runtime_error("Jensen-Shannon distance: element in p is < 0");
	}
	if (find_if(q.begin(), q.end(), checkNeg) != q.end()) {
		throw runtime_error("Jensen-Shannon distance: element in q is < 0");
	}
	double tol = 0.001;
	size_t expected = p.shape()[0];
	if (p.shape().size() == 1) {
		expected = 1;
	}
	auto s1 = xt::sum(p)[0];
	if (fabs(s1 - expected) > tol) {
		throw runtime_error("Jensen-Shannon distance: values in p do not sum to 1");
	}
	auto s2 = xt::sum(q)[0];
	if (fabs(s2 - expected) > tol) {
		throw runtime_error("Jensen-Shannon distance: values in q do not sum to 1");
	}
}

double jensen_shannon_distance(
	const xt::xarray<double>& p, 
	const xt::xarray<double>& q
) {
	jensen_shannon_distance_argument_check(p, q);
	xt::xarray<double> m = (p + q) / 2;
	double divergence = (
		_kl_divergence_double_ec(p, m, false) + 
		_kl_divergence_double_ec(q, m, false)
	) / 2;
	return sqrt(max(divergence, (double) 0));
}

xt::xarray<double> jensen_shannon_distance(
	const xt::xarray<double>& p, 
	const xt::xarray<double>& q,
	const size_t axis
) {
	jensen_shannon_distance_argument_check(p, q);
	xt::xarray<double> m = (p + q) / 2;
	return xt::eval(
		xt::sqrt(
			xt::clip(
				(
					_kl_divergence_array_ec(p, m, axis, false) + 
					_kl_divergence_array_ec(q, m, axis, false)
				) / 2,	
				(double) 0.0,
				numeric_limits<double>::max()
			)
		)
	);
}
