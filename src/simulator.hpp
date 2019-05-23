#pragma once
#include "classical_density_predictor.hpp"
#include "convolution.hpp"
#include "poisson_solver.hpp"
#include "quantum_density_predictor.hpp"
#include "schrodinger_solver.hpp"

#include <es_fe/geometry.hpp>
#include <es_fe/mesh/mesh1.hpp>
#include <es_la/dense.hpp>
#include <es_la/io.hpp>
#include <es_util/numeric.hpp>
#include <es_util/phys.hpp>

#include <mkl.h>
#include <mkl_dfti.h>

#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>

#include <fstream>

// Converts FWHM to sigma for the normal distribution
double from_gauss_fwhm(double width)
{
	return width / (2 * std::sqrt(2 * std::log(2)));
}

class Simulator
{
public:
	Simulator()
	{}

	void run(int /* argc */, const char** /* argv */)
	{
		using namespace es_util::au::literals;

		///////////////////////////////////////////////////////////////////////
		//* Parameters */

		Params params;
		params.length = 200_nm;
		params.lattice_temp = 10_kelvin;
		params.effective_mass = .2;
		params.eps = 10;
		// params.dopant_conc = 2e20_per_cm3;
		// params.ec_surf = 0.7_evolt;
		params.dopant_conc = 7.7e19_per_cm3;
		params.ec_surf = .7_evolt;
		//params.ec_surf = 0_evolt;

		es_fe::Linear_grid grid;
		grid.add_tick(0_nm);
		grid.add_tick(1000, params.length);

		///////////////////////////////////////////////////////////////////////
		//* Quasi-classical solution */

		const auto x_grid = grid.grid();
		es_fe::Mesh1 mesh(x_grid);

		Poisson_solver<Classical_density_predictor> cl_solver(mesh, params);
		cl_solver.init();
		cl_solver.solve();
		cl_solver.write("p0.mat");

		///////////////////////////////////////////////////////////////////////
		//* Quantum solution */

		Poisson_solver<Quantum_density_predictor> q_solver(mesh, params);
		q_solver.init();
		q_solver.set_init_guess(cl_solver.solution());

		Schrodinger_solver schrod_solver(mesh, params, q_solver.solution_view());
		q_solver.density_predictor().set_schrodinger_view(schrod_solver.solution_view());

		schrod_solver.init();
		schrod_solver.solve();
		schrod_solver.write("q.mat");

		//return;

		q_solver.solve();
		q_solver.write("p1.mat");

		for (int i = 0; i < 5; ++i)
		{
			schrod_solver.solve();
			q_solver.solve();
		}

		q_solver.write("p2.mat");

		[[maybe_unused]] auto psi = schrod_solver.solution_view();

		//return;

		//////////////////////////////////////////////////////////////////////
		//* Calculation of ARPES spetrum */

		std::size_t ne = 1000;
		const auto e_min = -.8_evolt;
		const auto e_max = .6_evolt;

		std::size_t nk = 1000;
		const auto kx_min = -.3 / 1_ang;
		const auto kx_max = .3 / 1_ang;

		// const auto de_disorder = 0.125_evolt;
		// const auto de_app = from_gauss_fwhm(.265_evolt);
		// const auto dkx_app = from_gauss_fwhm(0.023 / 1_ang);

		const auto de_disorder = 0.075_evolt;
		const auto de_app = from_gauss_fwhm(.05_evolt);
		const auto dkx_app = from_gauss_fwhm(0.01 / 1_ang);

		const auto lambda = 5_nm;

		const auto es = es_util::Linear_grid<double>::from_min_max(e_min, e_max, ne);
		const auto ks = es_util::Linear_grid<double>::from_min_max(kx_min, kx_max, nk);

		const auto fermi = charge_neutral_fermi_level(params);

		// auto poisson_solution = q_solver.solution_view();
		//auto poisson_solution = cl_solver.solution_view();

		//////////////////////////////////////////////////////////////////////

		// const auto N = grid.size();
		// const auto Nk = N / 2 + 1;

		// DFTI_DESCRIPTOR_HANDLE fft;
		// auto st = DftiCreateDescriptor(&fft, DFTI_DOUBLE, DFTI_REAL, 1, N);
		// st = DftiSetValue(fft, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
		// st = DftiSetValue(fft, DFTI_NUMBER_OF_TRANSFORMS, psi.size());
		// st = DftiSetValue(fft, DFTI_INPUT_DISTANCE, N);
		// st = DftiSetValue(fft, DFTI_OUTPUT_DISTANCE, Nk);
		// st = DftiSetValue(fft, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX);
		// st = DftiCommitDescriptor(fft);

		// es_la::Matrix_xd psi_n(N, psi.size(), 0);
		// for (std::size_t j = 0; j < psi.size(); ++j)
		// 	for (std::size_t i = 0; i < grid.size() - 2; ++i)
		// 		psi_n(i, j) = psi(i, j) * std::exp(-zs(i) / lambda);

		// es_la::Matrix_x<std::complex<double>> psi_k(Nk, psi.size());

		// st = DftiComputeForward(fft, psi_n.data(), psi_k.data());

		// es_la::Matrix_x<double> psi_k_sq(Nk, psi.size());
		// for (std::size_t j = 0; j < psi.size(); ++j)
		// 	for (std::size_t i = 0; i < Nk; ++i)
		// 		psi_k_sq(i, j) = std::norm(psi_k(i, j));

		// {
		// 	es_la::Matfile_writer mf("qf.mat");
		// 	mf.write("psi_z", psi_n);
		// 	mf.write("psi_k", psi_k_sq);
		// }

		// DftiFreeDescriptor(&fft);

		//		return;

		//////////////////////////////////////////////////////////////////////

		es_la::Vector_xd int_exp_psi(psi.size());
		for (std::size_t i = 0; i < psi.size(); ++i)
		{
			const auto zs = [&](auto v) { return mesh.vertex(v).x(); };
			const auto fn = [&](auto v) { return std::exp(-zs(v) / lambda) * psi.at(v, i); };
			int_exp_psi[i] = es_util::trapez_int(mesh.n_vertices(), zs, fn);
		}

		es_la::Matrix_xd arp(ne, nk, 0);
		for (std::size_t ie = 0; ie < ne; ++ie)
		{
			const auto fd = es_util::fermi((es[ie] - fermi) / params.lattice_temp);
			for (std::size_t ik = 0; ik < nk; ++ik)
			{
				const auto k_sq_over_2m = es_util::sq(ks[ik]) / (2 * params.effective_mass);
				for (std::size_t i = 0; i < psi.size(); ++i)
				{
					const auto e = es[ie] - (psi[i] + k_sq_over_2m);
					auto lor = 1 / (1 + es_util::sq(e / de_disorder));
					arp(ie, ik) += lor * fd * es_util::sq(int_exp_psi[i]);
				}
			}
		}

		es_la::Matfile_writer mat("arp.mat");

		mat.write("e_min", es_util::au::to_evolt(e_min));
		mat.write("e_max", es_util::au::to_evolt(e_max));
		mat.write("kx_min", 1 / es_util::au::to_ang(1 / kx_min));
		mat.write("kx_max", 1 / es_util::au::to_ang(1 / kx_max));

		mat.write("de_disorder", es_util::au::to_evolt(de_disorder));
		mat.write("de_app", es_util::au::to_evolt(de_app));
		mat.write("dkx_app", 1 / es_util::au::to_ang(1 / dkx_app));

		mat.write("lambda", es_util::au::to_nm(lambda));
		mat.write("arp0", arp);

		gauss_rows_convolution(arp, dkx_app / (kx_max - kx_min));
		gauss_cols_convolution(arp, de_app / (e_max - e_min));

		mat.write("arp", arp);

		es_la::write_gnuplot_binary(
			"arp.dat", arp, [&ks](auto i) { return 1 / es_util::au::to_ang(1 / ks[i]); },
			[&es](auto i) { return es_util::au::to_evolt(es[i]); });

		// 		es_la::Matrix_xd arp_at_z(ne, nk);
		// 		auto arp = es_util::trapez_int(
		// 			grid.size(), zs,
		// 			[&](auto iz) {
		// 				arp_at_z = 0;

		// 				auto dump = std::exp(-zs(iz) / lambda);
		// 				if (dump < 1e-5)
		// 					return 0. * arp_at_z;

		// 				const auto phi = poisson_solution[static_cast<es_fe::Vertex_index>(iz)];

		// 				for (std::size_t ie = 0; ie < ne; ++ie)
		// 				{
		// 					const auto e = es[ie];

		// 					for (std::size_t ik = 0; ik < nk; ++ik)
		// 					{
		// //#define QM
		// #ifdef QM
		// 						for (std::size_t iq = 0; iq < psi.size(); ++iq)
		// 						{
		// 							const auto en = psi[iq];

		// 							const auto k = ks[ik];
		// 							const auto ek = es_util::sq(k) / (2 * params.effective_mass) - phi;

		// 							//const auto dir = es_util::fermi((e - fermi) / params.lattice_temp);
		// 							const auto dir = es_util::fermi((e - fermi) / params.lattice_temp);
		// 							//const auto dir = es_util::fd_int_half(-(e - fermi) / params.lattice_temp);

		// 							//const auto ee = e - ek;
		// 							const auto ee = e - en - es_util::sq(k) / (2 * params.effective_mass);
		// 							auto lor = 1 / (1 + es_util::sq(ee / delta_e));
		// 							auto psi_sq = 1;//psi(iz, iq) * psi(iz, iq);

		// 							arp_at_z(ie, ik) += lor * dir * psi_sq;
		// 						}
		// #else
		// 						const auto k = ks[ik];
		// 						const auto ek = es_util::sq(k) / (2 * params.effective_mass) - phi;

		// 						// const auto dir = es_util::fermi((e - fermi) / params.lattice_temp);
		// 						const auto dir = es_util::fd_int_minus_half(-(e - fermi) / params.lattice_temp);

		// 						const auto ee = e - ek;
		// 						//const auto ee = e - en - es_util::sq(k) / (2 * params.effective_mass);
		// 						auto lor = 1 / (1 + es_util::sq(ee / delta_e));
		// 						arp_at_z(ie, ik) += lor * dir;
		// #endif
		// 					}
		// 				}

		// 				std::cout << dump << '\n';
		// 				return dump * arp_at_z;
		// 			},
		// 			es_la::Matrix_xd(ne, nk, 0));

		// 		es_la::Matfile_writer mat("arp.mat");
		// 		mat.write("arp", arp);
	}
};
