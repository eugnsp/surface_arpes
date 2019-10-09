#pragma once
#include <esl/dense.hpp>

#include <mkl_dfti.h>
#include <mkl_types.h>

#include <complex>
#include <stdexcept>

#define CALL_MKL_DFTI(fn, ...)                                                                                         \
	do                                                                                                                 \
	{                                                                                                                  \
		const auto status = fn(__VA_ARGS__);                                                                           \
		if (status != 0)                                                                                               \
			throw std::runtime_error(::DftiErrorMessage(status));                                                      \
	} while (0)

inline esl::Matrix_x<std::complex<double>> fft_1d_cols_real_to_half_complex(const esl::Matrix_xd& in)
{
	const auto fft_size = in.rows() / 2 + 1;
	esl::Matrix_x<std::complex<double>> out(fft_size, in.cols());

	::DFTI_DESCRIPTOR_HANDLE handle;
	CALL_MKL_DFTI(DftiCreateDescriptor, &handle, DFTI_DOUBLE, DFTI_REAL, 1, in.rows());
	CALL_MKL_DFTI(::DftiSetValue, handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
	CALL_MKL_DFTI(::DftiSetValue, handle, DFTI_NUMBER_OF_TRANSFORMS, in.cols());
	CALL_MKL_DFTI(::DftiSetValue, handle, DFTI_INPUT_DISTANCE, in.rows());
	CALL_MKL_DFTI(::DftiSetValue, handle, DFTI_OUTPUT_DISTANCE, fft_size);
	CALL_MKL_DFTI(::DftiSetValue, handle, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX);
	CALL_MKL_DFTI(::DftiCommitDescriptor, handle);
	CALL_MKL_DFTI(::DftiComputeForward, handle, const_cast<double*>(in.data()), out.data());
	CALL_MKL_DFTI(::DftiFreeDescriptor, &handle);

	return out;
}
