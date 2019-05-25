#pragma once

#include <cstddef>

struct Params
{
	double length;
	std::size_t nx;
	double eps;
	double temp;
	double m_eff;
	double dopant_conc;
	double ec_surf;
	double de_disorder;

	std::size_t n_max_iters;
	double stop_ec_sup_norm;

	double de_inst;
	double dkx_inst;
	double mfp;

	double e_min;
	double e_max;
	std::size_t ne;

	double kx_max;
	std::size_t nkx;

	double kz_max;
};
