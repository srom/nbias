/*
Implementation of probability's product rule (for independent variables) 
and marginalization.
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

double product_rule_log(const xt::xarray<double>& log_p) {
	return xt::sum(log_p)[0];
}

xt::xarray<double> product_rule_log(const xt::xarray<double>& log_p, const size_t axis) {
	return xt::eval(xt::sum(log_p, std::move(axis)));
}

double product_rule_normalized(const xt::xarray<double>& p) {
	auto n_elements = p.shape()[0];
	return xt::exp(xt::sum(xt::log(p)) / n_elements)[0];
}

xt::xarray<double> product_rule_normalized(
	const xt::xarray<double>& p, 
	const size_t axis
) {
	auto n_elements = p.shape()[axis];
	return xt::eval(xt::exp(xt::sum(xt::log(p), std::move(axis)) / n_elements));
}

double log_sum_exp(const xt::xarray<double>& log_p) {
	double max_val = xt::amax(log_p)[0];
	return max_val + xt::log(xt::sum(xt::exp(log_p - max_val)))[0];
}

xt::xarray<double> log_sum_exp(const xt::xarray<double>& log_p, const size_t axis) {
	size_t ax = axis;
	xt::xarray<double> max_per_ax = xt::amax(log_p, std::move(ax));
	xt::xarray<double> max_val = xt::zeros<double>(log_p.shape());
	for (int i = 0; i < log_p.shape()[0] ; ++i) {
		for (int j = 0; j < log_p.shape()[1] ; ++j) {
			max_val(i, j) = max_per_ax[i];
		}
	}
	return xt::eval(
		max_per_ax + xt::log(xt::sum(xt::exp(log_p - max_val), std::move(axis)))
	);
}

double marginalization(
	const xt::xarray<double>& prior, 
	const xt::xarray<double>& likelihood
) {
	return xt::sum(xt::exp(xt::log(prior) + xt::log(likelihood)))[0];
}

xt::xarray<double> marginalization(
	const xt::xarray<double>& prior, 
	const xt::xarray<double>& likelihood, 
	const size_t axis
) {
	return xt::eval(
		xt::sum(
			xt::exp(xt::log(prior) + xt::log(likelihood)), 
			std::move(axis)
		)
	);
}

double marginalization_log(
	const xt::xarray<double>& log_prior, 
	const xt::xarray<double>& log_likelihood
) {
	return log_sum_exp(log_prior + log_likelihood);
}

xt::xarray<double> marginalization_log(
	const xt::xarray<double>& log_prior, 
	const xt::xarray<double>& log_likelihood, 
	const size_t axis
) {
	return log_sum_exp(log_prior + log_likelihood, axis);
}
