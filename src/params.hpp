#pragma once
#include <es_util/numeric.hpp>

#include <cmath>

struct Params
{
	double length;

	double lattice_temp;
	double effective_mass;
	double dopant_conc;
	double ec_surf;
};

inline double effective_2d_dos(const Params& p)
{
	return p.effective_mass * p.lattice_temp / es_util::math::pi;
}

inline double effective_3d_dos(const Params& p)
{
	return 2 * std::pow(p.effective_mass * p.lattice_temp / es_util::math::two_pi, 1.5);
}

inline double charge_neutral_fermi_level(const Params& p)
{
	const auto e_dos = effective_3d_dos(p);
	return p.lattice_temp * es_util::inverse_fd_int_half(p.dopant_conc / e_dos);
}
