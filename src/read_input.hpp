#pragma once
#include "params.hpp"
#include "tools.hpp"

#include <es_util/phys.hpp>

#include <istream>
#include <sstream>
#include <stdexcept>
#include <string>

[[nodiscard]] Params read_input(std::istream& input)
{
	Params p;

	const unsigned int n_params = 17;
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
				p.nx = static_cast<std::size_t>(d);
				break;

			case 2:
				p.temp = es_util::au::from_kelvin(d);
				break;

			case 3:
				p.m_eff = d;
				break;

			case 4:
				p.eps = d;
				break;

			case 5:
				p.dopant_conc = es_util::au::from_per_cm3(d);
				break;

			case 6:
				p.ec_surf = es_util::au::from_evolt(d);
				break;

			case 7:
				p.de_disorder = from_lorentz_fwhm(es_util::au::from_evolt(d));
				break;

			case 8:
				p.de_inst = from_gauss_fwhm(es_util::au::from_evolt(d));
				break;

			case 9:
				p.dkx_inst = from_gauss_fwhm(es_util::au::from_per_ang(d));
				break;

			case 10:
				p.mfp = es_util::au::from_nm(d);
				break;

			case 11:
				p.e_min = es_util::au::from_evolt(d);
				break;

			case 12:
				p.e_max = es_util::au::from_evolt(d);
				break;

			case 13:
				p.ne = static_cast<std::size_t>(d);
				break;

			case 14:
				p.kx_max = es_util::au::from_per_ang(d);
				break;

			case 15:
				p.nkx = static_cast<std::size_t>(d);
				break;

			case 16:
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
