/*
Implementation of KL divergence and Jensen-Shannon distance.
*/
#include <iostream>
#include <algorithm>
#include <cmath>
#include <xtensor/xarray.hpp>
#include <xtensor/xmath.hpp>

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
	size_t div = 1;
	if (p.shape().size() > 1) {
		div = p.shape()[0];
	}
	auto s1 = xt::sum(p)[0];
	if (fabs((s1 / (double) div) - 1) > tol) {
		throw runtime_error("KL divergence: values in p do not sum to 1");
	}
	auto s2 = xt::sum(q)[0];
	if (fabs((s2 / (double) div) - 1) > tol) {
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
	return xt::sum(p * (xt::nan_to_num(xt::log2(p)) - xt::log2(q)))[0];
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
	return xt::eval(xt::sum(p * (xt::nan_to_num(xt::log2(p)) - xt::log2(q)), move(axis)));
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
	size_t div = 1;
	if (p.shape().size() > 1) {
		div = p.shape()[0];
	}
	auto s1 = xt::sum(p)[0];
	if (fabs((s1 / (double) div) - 1) > tol) {
		throw runtime_error("Jensen-Shannon distance: values in p do not sum to 1");
	}
	auto s2 = xt::sum(q)[0];
	if (fabs((s2 / (double) div) - 1) > tol) {
		throw runtime_error("Jensen-Shannon distance: values in q do not sum to 1");
	}
}

xt::xarray<double> duplicate_rows(const xt::xarray<double>& a, const int n_rows) {
	const int n_cols = a.size();
	xt::xarray<double> out = xt::zeros<double>({n_rows, n_cols});
	for (int i = 0; i < n_rows; ++i) {
		for (int j = 0; j < n_cols; ++j) {
			out(i, j) = a[j];
		}
	}
	return out;
}
