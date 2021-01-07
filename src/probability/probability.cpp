/*
Implementation of the probability's product rule and marginalization for independent variables.
*/
#include <iostream>
#include <algorithm>
#include <cmath>
#include <xtensor/xarray.hpp>
#include <xtensor/xmath.hpp>

using namespace std;

double product_rule(const xt::xarray<double> p) {
	auto n_elements = p.shape()[0];
	return xt::exp(xt::sum(xt::log(p)) / n_elements)[0];
}

xt::xarray<double> product_rule(const xt::xarray<double> p, const size_t axis) {
	auto n_elements = p.shape()[axis];
	return xt::eval(xt::exp(xt::sum(xt::log(p), move(axis)) / n_elements));
}

double marginalization(const xt::xarray<double> prior, const xt::xarray<double> likelihood) {
	return xt::sum(xt::exp(xt::log(prior) + xt::log(likelihood)))[0];
}

xt::xarray<double> marginalization(
	const xt::xarray<double> prior, 
	const xt::xarray<double> likelihood, 
	const size_t axis
) {
	return xt::eval(xt::sum(xt::exp(xt::log(prior) + xt::log(likelihood)), move(axis)));
}
