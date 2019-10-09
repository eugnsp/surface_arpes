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
#include <esl/dense.hpp>
#include <esl/io.hpp>
#include <esu/numeric.hpp>
#include <esu/phys.hpp>

#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>

class Simulator
{
public:
	Simulator()
	{}

	auto read_params(const char* file_name = nullptr)
	{
		if (!file_name)
			return read_input(std::cin);
		else
			return read_input(std::ifstream(file_name));
	}

	void run([[maybe_unused]] int argc, [[maybe_unused]] const char** argv)
	{
		///////////////////////////////////////////////////////////////////////
		//* Parameters */

#ifdef NDEBUG
		const auto p = (argc == 2) ? read_params(argv[1]) : read_params();
#else
		const auto p = read_params("input.txt");
#endif

		es_fe::Linear_grid grid;
		grid.add_tick(0);
		grid.add_tick(p.nx - 1, p.length);

		///////////////////////////////////////////////////////////////////////
		//* The Poisson-Schrodinger solver */

		const auto fermi = charge_neutral_fermi_level(p);

		std::cout << "----------------------------------------------\n"
				  << "Running the Poisson-Schrodinger solver\n"
				  << std::endl;

		std::cout << "System length: L = " << esu::au::to_nm(p.length) << " nm\n"
				  << "Temperature: T = " << esu::au::to_kelvin(p.temp) << " K\n"
				  << "Effective mass: m_eff = " << p.m_eff << "\n"
				  << "Static dielectric permittivity: eps = " << p.eps << "\n"
				  << "Dopant concentration: N_d = " << esu::au::to_per_cm3(p.dopant_conc) << " cm^-3\n"
				  << "Surface potential: Ec(0) = " << esu::au::to_evolt(p.ec_surf) << " eV\n"
				  << "Disorder strength: gamma_D = " << esu::au::to_evolt(p.gamma_disorder) << " eV\n"
				  << "Fermi level: F = " << esu::au::to_evolt(fermi) << " eV\n"
				  << "Number of potential mesh points: N_x = " << p.nx << "\n"
				  << "Poisson-Schrodinger solver maximum number of iterations: N_max = " << p.n_max_iters << "\n"
				  << "Poisson-Schrodinger solver stopping criterion: |phi - phi_prev|_s < "
				  << esu::au::to_evolt(p.stop_ec_sup_norm) << " eV\n"
				  << std::endl;

		///////////////////////////////////////////////////////////////////////
		//* Quasi-classical solution */

		es_fe::Mesh1 mesh(grid.grid());

		Poisson_solver<Classical_density_predictor> cl_solver(mesh, p);
		cl_solver.init();
		cl_solver.solve();
		cl_solver.write_mat("poisson_cl.mat");

		///////////////////////////////////////////////////////////////////////
		//* Quantum solution */

		Poisson_solver<Quantum_density_predictor> q_solver(mesh, p);
		q_solver.init();
		q_solver.set_init_guess(cl_solver.solution());

		Schrodinger_solver schrod_solver(mesh, p, q_solver.solution_view());
		q_solver.density_predictor().set_schrodinger_view(schrod_solver.solution_view());

		schrod_solver.init();
		for (std::size_t i = 0;;)
		{
			schrod_solver.solve();
			q_solver.solve();

			auto sup_norm = q_solver.density_predictor().potential_change_sup_norm();
			std::cout << i + 1 << ". Potential change |phi - phi_prev|_s = " << esu::au::to_evolt(sup_norm) << " eV"
					  << std::endl;

			if (sup_norm < p.stop_ec_sup_norm)
				break;
			if (++i >= p.n_max_iters)
				throw std::runtime_error("No convergence in the Poisson-Schrodinger solver");
		}
		q_solver.write_mat("poisson_q.mat");
		schrod_solver.write("schrod.mat");
		//q_solver.write_la("poisson_q.la");
		// q_solver.read_la("poisson_q.la");
		// schrod_solver.solve();

		//////////////////////////////////////////////////////////////////////
		//* Calculation of ARPES spetrum */

		std::cout << "----------------------------------------------\n"
				  << "Running ARPES spectrum calculator\n"
				  << std::endl;

		const auto psi = schrod_solver.solution_view();
		esl::Matrix_xd psi_z(*mesh.n_vertices(), psi.size());
		for (auto& v : mesh.vertices())
		{
			const auto exp = std::exp(-v.vertex().x() / p.mfp);
			for (std::size_t j = 0; j < psi.size(); ++j)
				psi_z(**v, j) = exp * psi.at(*v, j);
		}

		const auto psi_k = fft_1d_cols_real_to_half_complex(psi_z);

		const auto dkz = esu::math::two_pi / p.length;
		const auto nkz = std::min(static_cast<std::size_t>(std::ceil(p.kz_max / dkz)), psi_k.rows());

		const auto ikx_max = p.nkx - 1;
		const auto nkx_f = 2 * ikx_max + 1;

		const auto ikz_max = nkz - 1;
		const auto nkz_f = 2 * ikz_max + 1;

		const auto kxs = esu::Linear_grid<double>::from_min_max(0, p.kx_max, p.nkx);
		const auto kzs = esu::Linear_grid<double>::from_min_step(0, dkz, nkz);
		const auto es = esu::Linear_grid<double>::from_min_max(p.e_min, p.e_max, p.ne);

		//////////////////////////////////////////////////////////////////////

		std::cout << "Instrumental energy broadening: sigma_E = " << esu::au::to_evolt(p.sigma_e_inst) << " eV\n"
				  << "Instrumental Kx broadening: sigma_Kx = " << esu::au::to_per_ang(p.sigma_kx_inst)
				  << " Ang^-1\n"
				  << "Electron's mean free path: lambda = " << esu::au::to_nm(p.mfp) << " nm\n"
				  << "Minimum energy: E_min = " << esu::au::to_evolt(p.e_min) << " eV\n"
				  << "Maximum energy: E_max = " << esu::au::to_evolt(p.e_max) << " eV\n"
				  << "Energy resolution: N_E = " << p.ne << "\n"
				  << "Energy resolution: delta_E = " << esu::au::to_evolt((p.e_max - p.e_min) / p.ne) << " eV\n"
				  << "Maximum Kx: Kx_max = " << esu::au::to_per_ang(p.kx_max) << " Ang^-1\n"
				  << "Kx mesh size: N_Kx = " << nkx_f << "\n"
				  << "Kx resolution: delta_Kx = " << esu::au::to_per_ang(p.kx_max / p.nkx) << " Ang^-1\n"
				  << "Maximum Kz: Kz_max = " << esu::au::to_per_ang(kzs.back()) << " Ang^-1\n"
				  << "Kz mesh size: N_Kz = " << nkz_f << "\n"
				  << "Kz resolution: delta_Kz = " << esu::au::to_per_ang(dkz) << " Ang^-1\n"
				  << std::endl;

		//////////////////////////////////////////////////////////////////////

		esl::Matrix_xd arp_kx_e(nkx_f, p.ne, 0);
		esl::Matrix_xd arp_kx_kz(nkx_f, nkz_f, 0);
		esl::Matrix_xd arp_e_kz(p.ne, nkz_f, 0);

		const auto ie_fermi = static_cast<std::size_t>(std::round((fermi - es[0]) / (es[1] - es[0])));
		if (ie_fermi >= p.ne)
			throw std::runtime_error("Fermi level is outside the specified energy range");

		esl::Matrix_xd arp_n(nkx_f, p.ne);
		for (std::size_t ip = 0; ip < psi.size(); ++ip)
		{
			for (std::size_t ie = 0; ie < p.ne; ++ie)
			{
				const auto f_fd = esu::fermi((es[ie] - fermi) / p.temp);
				for (std::size_t ikx = 0; ikx < p.nkx; ++ikx)
				{
					const auto k_sq_over_2m = esu::sq(kxs[ikx]) / (2 * p.m_eff);
					const auto e = es[ie] - (psi[ip] + k_sq_over_2m);
					const auto f_d = 1 / (1 + esu::sq(e / p.gamma_disorder));
					arp_n(ikx_max + ikx, ie) = arp_n(ikx_max - ikx, ie) = f_d * f_fd;
				}
			}

			gauss_cols_convolution(arp_n, p.sigma_kx_inst / (2 * p.kx_max));
			gauss_rows_convolution(arp_n, p.sigma_e_inst / (p.e_max - p.e_min));

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

		esl::Matfile_writer mat("arpes.mat");

		mat.write("e_min", esu::au::to_evolt(p.e_min));
		mat.write("e_max", esu::au::to_evolt(p.e_max));
		mat.write("kx_max", esu::au::to_per_ang(p.kx_max));
		mat.write("kz_max", esu::au::to_per_ang(kzs.back()));

		mat.write("gamma_disorder", esu::au::to_evolt(p.gamma_disorder));
		mat.write("sigma_e_inst", esu::au::to_evolt(p.sigma_e_inst));
		mat.write("sigma_kx_inst", esu::au::to_per_ang(p.sigma_kx_inst));

		mat.write("mfp", esu::au::to_nm(p.mfp));

		mat.write("arpes_kx_e", arp_kx_e);
		mat.write("arpes_kx_kz", arp_kx_kz);
		mat.write("arpes_e_kz", arp_e_kz);

		const auto kxs_title = [&kxs, ikx_max](auto i) {
			const auto kx = (i <= ikx_max) ? -kxs[ikx_max - i] : kxs[i - ikx_max];
			return esu::au::to_per_ang(kx);
		};

		const auto kzs_title = [&kzs, ikz_max](auto i) {
			const auto kz = (i <= ikz_max) ? -kzs[ikz_max - i] : kzs[i - ikz_max];
			return esu::au::to_per_ang(kz);
		};

		const auto es_title = [&es](auto i) { return esu::au::to_evolt(es[i]); };

		esl::write_gnuplot_binary("arpes_kx_e.dat", arp_kx_e, kxs_title, es_title);
		esl::write_gnuplot_binary("arpes_kx_kz.dat", arp_kx_kz, kxs_title, kzs_title);
		esl::write_gnuplot_binary("arpes_e_kz.dat", arp_e_kz, es_title, kzs_title);
	}
};
