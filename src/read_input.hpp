#pragma once
#include "params.hpp"
#include "tools.hpp"

#include <es_util/phys.hpp>

#include <istream>
#include <sstream>
#include <stdexcept>
#include <string>

[[nodiscard]] inline Params read_input(std::istream& input)
{
	Params p;

	const unsigned int n_params = 19;
	unsigned int i = 0;

	std::string str;
	while (std::getline(input, str) && i < n_params)
	{
		std::istringstream ss(str);

		double d;
		if (ss >> d)
		{
			switch (i)
			{
			case 0:
				p.length = es_util::au::from_nm(d);
				break;

			case 1:
				p.temp = es_util::au::from_kelvin(d);
				break;

			case 2:
				p.m_eff = d;
				break;

			case 3:
				p.eps = d;
				break;

			case 4:
				p.dopant_conc = es_util::au::from_per_cm3(d);
				break;

			case 5:
				p.ec_surf = es_util::au::from_evolt(d);
				break;

			case 6:
				p.gamma_disorder = from_lorentz_fwhm(es_util::au::from_evolt(d));
				break;

			case 7:
				p.nx = static_cast<std::size_t>(d);
				break;

			case 8:
				p.n_max_iters = static_cast<std::size_t>(d);
				break;

			case 9:
				p.stop_ec_sup_norm = es_util::au::from_evolt(d);
				break;

			case 10:
				p.sigma_e_inst = from_gauss_fwhm(es_util::au::from_evolt(d));
				break;

			case 11:
				p.sigma_kx_inst = from_gauss_fwhm(es_util::au::from_per_ang(d));
				break;

			case 12:
				p.mfp = es_util::au::from_nm(d);
				break;

			case 13:
				p.e_min = es_util::au::from_evolt(d);
				break;

			case 14:
				p.e_max = es_util::au::from_evolt(d);
				break;

			case 15:
				p.ne = static_cast<std::size_t>(d);
				break;

			case 16:
				p.kx_max = es_util::au::from_per_ang(d);
				break;

			case 17:
				p.nkx = static_cast<std::size_t>(d);
				break;

			case 18:
				p.kz_max = es_util::au::from_per_ang(d);
				break;
			}

			++i;
		}
	}

	if (i < n_params)
		throw std::runtime_error("Bad configuration file");

	return p;
}
