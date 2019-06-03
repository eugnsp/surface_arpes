#pragma once
#include <es_la/dense.hpp>

#include <es_util/numeric.hpp>

#include <cmath>
#include <cstddef>

inline es_la::Vector_xd gauss_distribution(const std::size_t n, const double alpha)
{
	es_la::Vector_xd vec(n);

	const auto one_over_sigma = 1 / (alpha * n);
	const auto pre = one_over_sigma / es_util::math::sqrt_two_pi;
	for (std::size_t i = 0; i < n; ++i)
	{
		const auto z = i * one_over_sigma - .5 / alpha;
		vec[i] = pre * std::exp(-.5 * z * z);
	}

	return vec;
}

inline void gauss_cols_convolution(es_la::Matrix_xd& y, const double alpha)
{
	const auto x = gauss_distribution(y.rows(), alpha);
	es_la::cols_convolution(x, y, y.rows() / 2);
}

inline void gauss_rows_convolution(es_la::Matrix_xd& y, const double alpha)
{
	const auto x = gauss_distribution(y.cols(), alpha);
	es_la::rows_convolution(x, y, y.cols() / 2);
}
