#pragma once
#include "params.hpp"

#include <es_util/numeric.hpp>

#include <cmath>

inline double effective_2d_dos(const Params& p)
{
	return p.m_eff * p.temp / es_util::math::pi;
}

inline double effective_3d_dos(const Params& p)
{
	return 2 * std::pow(p.m_eff * p.temp / es_util::math::two_pi, 1.5);
}

inline double charge_neutral_fermi_level(const Params& p)
{
	const auto e_dos = effective_3d_dos(p);
	return p.temp * es_util::inverse_fd_int_half(p.dopant_conc / e_dos);
}

// Converts FWHM to sigma for the normal distribution
inline double from_gauss_fwhm(double width)
{
	return width / (2 * std::sqrt(2 * std::log(2)));
}

// Converts FWHM to gamma for the Lorentz distribution
inline double from_lorentz_fwhm(double width)
{
	return width / 2;
}
