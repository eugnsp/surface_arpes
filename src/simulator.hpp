#pragma once
#include "classical_density_predictor.hpp"
#include "poisson_solver.hpp"
#include "quantum_density_predictor.hpp"
#include "schrodinger_solver.hpp"

#include <es_fe/geometry.hpp>
#include <es_fe/mesh/mesh1.hpp>
#include <es_la/dense.hpp>
#include <es_la/io/matfile_writer.hpp>
#include <es_util/numeric.hpp>
#include <es_util/phys.hpp>

#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>

class Simulator
{
public:
	Simulator()
	{}

	void run(int /* argc */, const char** /* argv */)
	{
		using namespace es_util::au::literals;

		{
			es_la::Matrix_d<5, 5> m1;
			es_la::Matrix_d<5, 5> m2;
			auto m = m1+ m2;
			auto mm = m.view(1, 2, 1, 2);
			//auto z = mm(1, 1);
		}

		Params params;
		params.length = 50_nm;
		params.lattice_temp = 10_kelvin;
		params.effective_mass = .2;
		// params.dopant_conc = 2e20_per_cm3;
		// params.ec_surf = 0.7_evolt;
		params.dopant_conc = 7.7e19_per_cm3;
		params.ec_surf = 0.7_evolt;

		es_fe::Linear_grid grid;
		grid.add_tick(0_nm);
		grid.add_tick(100, params.length / 2, .9);
		grid.add_tick(200, params.length, -.9);

		const auto x_grid = grid.grid();
		es_fe::Mesh1 mesh_(x_grid);

		Classical_density_predictor cdp(params);

		Poisson_solver solver1(mesh_, params, cdp);
		solver1.init();
		solver1.solve();
		solver1.write("p0.mat");

		Quantum_density_predictor qdp(params);
		Poisson_solver solver2(mesh_, params, qdp);

		solver2.init();
		solver2.solution_view().raw() = solver1.solution_view().raw();

		Schrodinger_solver ss(mesh_, params, solver2.solution_view());
		ss.init();
		ss.solve();
		ss.write("q.mat");

		for (std::size_t i = 0; i < ss.solution_view().size(); ++i)
		{
			std::cout << ss.solution_view()[i] << std::endl;
		}

		qdp.set_schrodinger_view(ss.solution_view());

		qdp.set_potential_view(solver2.solution_view());
		qdp.xxx();

		//return;

		solver2.solve();
		solver2.write("p1.mat");

		for (int i = 0; i < 5; ++i)
		{
			qdp.xxx();
			ss.solve();
			solver2.solve();
		}

		qdp.xxx();
		solver2.write("p2.mat");

		auto psi = ss.solution_view();

		return;

		std::size_t ne = 500;
		std::size_t nk = 500;

		auto lambda = 1_nm;

		auto es = es_util::Linear_grid<double>::from_min_max(-1.4_evolt + .33_evolt, .45_evolt, ne);
		auto ks = es_util::Linear_grid<double>::from_min_max(-.3 / 1_ang, .3 / 1_ang, nk);

		const auto delta_e = .125_evolt;

		const auto fermi = charge_neutral_fermi_level(params);

		const auto zs = [&x_grid](auto i) { return x_grid[i].x(); };

		auto poisson_solution = solver2.solution_view();

		es_la::Matrix_xd arp_at_z(ne, nk);
		auto arp = es_util::trapez_int(
			grid.size(), zs,
			[&](auto iz) {
				arp_at_z.zero();

				auto dump = std::exp(-zs(iz) / lambda);
				if (dump < 1e-5)
					return 0. * arp_at_z;

				const auto phi = poisson_solution[static_cast<es_fe::Vertex_index>(iz)];

				for (std::size_t ie = 0; ie < ne; ++ie)
				{
					const auto e = es[ie];
					const auto dir = es_util::fermi((e - fermi) / params.lattice_temp);

					for (std::size_t ik = 0; ik < nk; ++ik)
					{
						// for (std::size_t iq = 0; iq < psi.size(); ++iq)
						// {
						// const auto en = psi[iq];

						const auto k = ks[ik];
						const auto ee = e + phi - es_util::sq(k) / (2 * params.effective_mass);
						//const auto ee = e - en - es_util::sq(k) / (2 * params.effective_mass);
						auto lor = 1 / (1 + es_util::sq(ee / delta_e));
						// auto psi_sq = psi(iz, iq) * psi(iz, iq);
						auto psi_sq = 1;

						arp_at_z(ie, ik) += lor * dir * psi_sq;
						// }
					}
				}

				std::cout << dump << '\n';
				return dump * arp_at_z;
			},
			es_la::Matrix_xd(ne, nk, 0));

		es_la::Matfile_writer mat("arp.mat");
		mat.write("arp", arp);

		for (std::size_t i = 0; i < psi.size(); ++i)
		{
			std::cout << psi[i] << '\n';
		}
	}
};
