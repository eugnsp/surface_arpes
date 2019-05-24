#pragma once
#include "classical_density_predictor.hpp"
#include "convolution.hpp"
#include "fft.hpp"
#include "params.hpp"
#include "poisson_solver.hpp"
#include "quantum_density_predictor.hpp"
#include "read_input.hpp"
#include "schrodinger_solver.hpp"
#include "tools.hpp"

#include <es_fe/geometry.hpp>
#include <es_fe/mesh/mesh1.hpp>
#include <es_la/dense.hpp>
#include <es_la/io.hpp>
#include <es_util/numeric.hpp>

#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>

class Simulator
{
public:
	Simulator()
	{}

	void run(int /* argc */, const char** /* argv */)
	{
		///////////////////////////////////////////////////////////////////////
		//* Parameters */

		const auto p = read_input(std::cin);
		// std::ifstream config("input.txt");
		// const auto p = read_input(config);

		es_fe::Linear_grid grid;
		grid.add_tick(0);
		grid.add_tick(p.nx - 1, p.length);

		///////////////////////////////////////////////////////////////////////
		//* Quasi-classical solution */

		es_fe::Mesh1 mesh(grid.grid());

		Poisson_solver<Classical_density_predictor> cl_solver(mesh, p);
		cl_solver.init();
		cl_solver.solve();
		cl_solver.write("poisson_cl.mat");

		///////////////////////////////////////////////////////////////////////
		//* Quantum solution */

		Poisson_solver<Quantum_density_predictor> q_solver(mesh, p);
		q_solver.init();
		q_solver.set_init_guess(cl_solver.solution());

		Schrodinger_solver schrod_solver(mesh, p, q_solver.solution_view());
		q_solver.density_predictor().set_schrodinger_view(schrod_solver.solution_view());

		schrod_solver.init();
		for (int i = 0; i < 5; ++i)
		{
			schrod_solver.solve();
			q_solver.solve();
		}
		q_solver.write("poisson_q.mat");

		//////////////////////////////////////////////////////////////////////
		//* Calculation of ARPES spetrum */

		const auto psi = schrod_solver.solution_view();
		es_la::Matrix_xd psi_z(*mesh.n_vertices(), psi.size());
		for (auto& v : mesh.vertices())
		{
			const auto exp = std::exp(-v.vertex().x() / p.mfp);
			for (std::size_t j = 0; j < psi.size(); ++j)
				psi_z(**v, j) = exp * psi.at(*v, j);
		}

		const auto psi_k = fft_1d_cols_real_to_half_complex(psi_z);

		const auto dkz = es_util::math::two_pi / p.length;
		const auto nkz = std::min(static_cast<std::size_t>(std::ceil(p.kz_max / dkz)), psi_k.rows());

		const auto ikx_max = p.nkx - 1;
		const auto nkx_f = 2 * ikx_max + 1;

		const auto ikz_max = nkz - 1;
		const auto nkz_f = 2 * ikz_max + 1;

		const auto kxs = es_util::Linear_grid<double>::from_min_max(0, p.kx_max, p.nkx);
		const auto kzs = es_util::Linear_grid<double>::from_min_step(0, dkz, nkz);
		const auto es = es_util::Linear_grid<double>::from_min_max(p.e_min, p.e_max, p.ne);

		es_la::Matrix_xd arp_kx_e(nkx_f, p.ne, 0);
		es_la::Matrix_xd arp_kx_kz(nkx_f, nkz_f, 0);
		es_la::Matrix_xd arp_e_kz(p.ne, nkz_f, 0);

		const auto fermi = charge_neutral_fermi_level(p);
		const auto ie_fermi = static_cast<std::size_t>(std::round((fermi - es[0]) / (es[1] - es[0])));
		if (ie_fermi >= p.ne)
			throw std::runtime_error("Fermi level is outside the specified energy range");

		es_la::Matrix_xd arp_n(nkx_f, p.ne);
		for (std::size_t ip = 0; ip < psi.size(); ++ip)
		{
			for (std::size_t ie = 0; ie < p.ne; ++ie)
			{
				const auto f_fd = es_util::fermi((es[ie] - fermi) / p.temp);
				for (std::size_t ikx = 0; ikx < p.nkx; ++ikx)
				{
					const auto k_sq_over_2m = es_util::sq(kxs[ikx]) / (2 * p.m_eff);
					const auto e = es[ie] - (psi[ip] + k_sq_over_2m);
					const auto f_d = 1 / (1 + es_util::sq(e / p.de_disorder));
					arp_n(ikx_max + ikx, ie) = arp_n(ikx_max - ikx, ie) = f_d * f_fd;
				}
			}

			gauss_cols_convolution(arp_n, p.dkx_inst / (2 * p.kx_max));
			gauss_rows_convolution(arp_n, p.de_inst / (p.e_max - p.e_min));

			const auto psi_k0_sq = std::norm(psi_k(0, ip));
			for (std::size_t ie = 0; ie < p.ne; ++ie)
				for (std::size_t ikx = 0; ikx < nkx_f; ++ikx)
					arp_kx_e(ikx, ie) += arp_n(ikx, ie) * psi_k0_sq;

			for (std::size_t ikz = 0; ikz < nkz; ++ikz)
			{
				const auto psi_k_sq = std::norm(psi_k(ikz, ip));
				for (std::size_t ikx = 0; ikx < nkx_f; ++ikx)
					arp_kx_kz(ikx, ikz_max + ikz) += arp_n(ikx, ie_fermi) * psi_k_sq;
				for (std::size_t ie = 0; ie < p.ne; ++ie)
					arp_e_kz(ie, ikz_max + ikz) += arp_n(ikx_max, ie) * psi_k_sq;
			}
		}

		// TODO : use flip views
		for (std::size_t ikz = 1; ikz < nkz; ++ikz)
		{
			arp_kx_kz.col_view(ikz_max - ikz) = arp_kx_kz.col_view(ikz_max + ikz);
			arp_e_kz.col_view(ikz_max - ikz) = arp_e_kz.col_view(ikz_max + ikz);
		}

		//////////////////////////////////////////////////////////////////////
		//* Export results */

		es_la::Matfile_writer mat("arpes.mat");

		mat.write("e_min", es_util::au::to_evolt(p.e_min));
		mat.write("e_max", es_util::au::to_evolt(p.e_max));
		mat.write("kx_max", es_util::au::to_per_ang(p.kx_max));
		mat.write("kz_max", es_util::au::to_per_ang(kzs.back()));

		mat.write("de_disorder", es_util::au::to_evolt(p.de_disorder));
		mat.write("de_inst", es_util::au::to_evolt(p.de_inst));
		mat.write("dkx_inst", es_util::au::to_per_ang(p.dkx_inst));

		mat.write("mfp", es_util::au::to_nm(p.mfp));

		mat.write("arpes_kx_e", arp_kx_e);
		mat.write("arpes_kx_kz", arp_kx_kz);
		mat.write("arpes_e_kz", arp_e_kz);

		const auto kxs_title = [&kxs, ikx_max](auto i) {
			const auto kx = (i <= ikx_max) ? -kxs[ikx_max - i] : kxs[i - ikx_max];
			return es_util::au::to_per_ang(kx);
		};

		const auto kzs_title = [&kzs, ikz_max](auto i) {
			const auto kz = (i <= ikz_max) ? -kzs[ikz_max - i] : kzs[i - ikz_max];
			return es_util::au::to_per_ang(kz);
		};

		const auto es_title = [&es](auto i) { return es_util::au::to_evolt(es[i]); };

		es_la::write_gnuplot_binary("arpes_kx_e.dat", arp_kx_e, kxs_title, es_title);
		es_la::write_gnuplot_binary("arpes_kx_kz.dat", arp_kx_kz, kxs_title, kzs_title);
		es_la::write_gnuplot_binary("arpes_e_kz.dat", arp_e_kz, es_title, kzs_title);
	}
};
