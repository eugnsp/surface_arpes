#pragma once
#include "params.hpp"
#include "poisson_solver_base.hpp"
#include "schrodinger_solver_base.hpp"

#include <es_la/dense.hpp>
#include <es_util/numeric.hpp>

#include <cstddef>
#include <memory>
#include <utility>

// Calculator of electron density using the Scrodinger equation and the perturbation theory
class Quantum_density_predictor
{
private:
	using Poisson_solution_view = Poisson_solver_base::Solution_view;
	using Poisson_solution = Poisson_solver_base::Solution;
	using Schrodinger_solution_view = Schrodinger_solver_base::Solution_view;

public:
	Quantum_density_predictor(const Params& params, Poisson_solution_view phi) :
		params_(params), phi_(phi), prev_phi_(phi)
	{
		effective_dos_ = effective_2d_dos(params_);
		fermi_level_ = charge_neutral_fermi_level(params_);
	}

	void before_solve()
	{
		prev_phi_ = phi_;
	}

	void set_schrodinger_view(Schrodinger_solution_view psi)
	{
		psi_ = std::make_unique<Schrodinger_solution_view>(psi);
	}

	template<class Element, class Quadr, class Dofs>
	auto get(const Dofs& dofs, const es_fe::Mesh1::Edge_view& edge) const
	{
		es_la::Vector<std::pair<double, double>, Quadr::size> density{};
		const auto phis = es_fe::at_quadr<Element, Quadr>(phi_, dofs);
		const auto prev_phis = es_fe::at_quadr<Element, Quadr>(prev_phi_, dofs);
		const auto dphis = (phis - prev_phis).eval();

		for (std::size_t is = 0; is < psi_->size(); ++is)
		{
			const auto energy = (*psi_)[is];
			const auto psis = es_fe::at_quadr<Element, Quadr>(*psi_, is, edge);

			for (std::size_t iq = 0; iq < Quadr::size; ++iq)
			{
				// Note: eigenvectors are normalized with respect to the scalar product
				// <f | g> = <f| B |g>, for the eigenproblem A |f> = lambda B |f>.

				const auto z = (-energy + dphis[iq] + fermi_level_) / params_.lattice_temp;
				const auto f = es_util::ln_one_p_exp(z);
				const auto fd = es_util::fermi(-z);
				const auto psi_sq = es_util::sq(psis[iq]);

				density[iq].first -= effective_dos_ * psi_sq * f;
				density[iq].second -= effective_dos_ / params_.lattice_temp * psi_sq * fd;
			}
		}

		for (std::size_t iq = 0; iq < Quadr::size; ++iq)
			density[iq].first += params_.dopant_conc;

		return density;
	}

private:
	const Params params_;
	double effective_dos_;
	double fermi_level_;

	const Poisson_solution_view phi_;
	Poisson_solution prev_phi_;

	std::unique_ptr<Schrodinger_solution_view> psi_;
};
