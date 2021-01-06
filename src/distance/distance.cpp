/*
Implementation of KL divergence and Jensen-Shannon distance.
*/
#include <iostream>
#include <algorithm>
#include <cmath>
#include <xtensor/xarray.hpp>
#include <xtensor/xmath.hpp>
#include "distance_util.cpp"

using namespace std;

double kl_divergence(
	const xt::xarray<double>& p, 
	const xt::xarray<double>& q
) {
	return _kl_divergence_double_ec(p, q, true);
}

xt::xarray<double> kl_divergence(
	const xt::xarray<double>& p, 
	const xt::xarray<double>& q,
	const size_t axis
) {
	return _kl_divergence_array_ec(p, q, axis, true);
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
