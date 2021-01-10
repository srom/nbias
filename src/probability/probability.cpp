/*
Implementation of the probability's product rule and marginalization for independent variables.
*/
#include <utility>
#include <xtensor/xarray.hpp>
#include <xtensor/xmath.hpp>

double product_rule(const xt::xarray<double>& p) {
	return xt::exp(xt::sum(xt::log(p)))[0];
}

xt::xarray<double> product_rule(const xt::xarray<double>& p, const size_t axis) {
	return xt::eval(xt::exp(xt::sum(xt::log(p), std::move(axis))));
}

double product_rule_normalized(const xt::xarray<double>& p) {
	auto n_elements = p.shape()[0];
	return xt::exp(xt::sum(xt::log(p)) / n_elements)[0];
}

xt::xarray<double> product_rule_normalized(const xt::xarray<double>& p, const size_t axis) {
	auto n_elements = p.shape()[axis];
	return xt::eval(xt::exp(xt::sum(xt::log(p), std::move(axis)) / n_elements));
}

double marginalization(const xt::xarray<double>& prior, const xt::xarray<double>& likelihood) {
	return xt::sum(xt::exp(xt::log(prior) + xt::log(likelihood)))[0];
}

xt::xarray<double> marginalization(
	const xt::xarray<double>& prior, 
	const xt::xarray<double>& likelihood, 
	const size_t axis
) {
	return xt::eval(xt::sum(xt::exp(xt::log(prior) + xt::log(likelihood)), std::move(axis)));
}
