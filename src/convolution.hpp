#pragma once
#include <esl/dense.hpp>

#include <esu/numeric.hpp>

#include <cmath>
#include <cstddef>

inline esl::Vector_xd gauss_distribution(const std::size_t n, const double alpha)
{
	esl::Vector_xd vec(n);

	const auto one_over_sigma = 1 / (alpha * n);
	const auto pre = one_over_sigma / esu::math::sqrt_two_pi;
	for (std::size_t i = 0; i < n; ++i)
	{
		const auto z = i * one_over_sigma - .5 / alpha;
		vec[i] = pre * std::exp(-.5 * z * z);
	}

	return vec;
}

inline void gauss_cols_convolution(esl::Matrix_xd& y, const double alpha)
{
	const auto x = gauss_distribution(y.rows(), alpha);
	esl::cols_convolution(x, y, y.rows() / 2);
}

inline void gauss_rows_convolution(esl::Matrix_xd& y, const double alpha)
{
	const auto x = gauss_distribution(y.cols(), alpha);
	esl::rows_convolution(x, y, y.cols() / 2);
}
