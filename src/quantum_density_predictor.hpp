#pragma once
#include "params.hpp"
#include "poisson_solver_base.hpp"
#include "schrodinger_solver_base.hpp"
#include "tools.hpp"

#include <esl/dense.hpp>
#include <esu/numeric.hpp>

#include <cstddef>
#include <optional>
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
		psi_.emplace(std::move(psi));
	}

	template<class Quadr, class Dofs>
	auto get(const Dofs& dofs, const esf::Mesh1::Edge_view& edge) const
	{
		std::pair<esl::Vector_d<Quadr::size>, esl::Vector_d<Quadr::size>> density;

		const auto phis = esf::at_quadr<Quadr>(phi_, dofs);
		const auto prev_phis = esf::at_quadr<Quadr>(prev_phi_, dofs);
		const auto dphis = (phis - prev_phis).eval();

		for (std::size_t is = 0; is < psi_->size(); ++is)
		{
			const auto energy = (*psi_)[is];
			const auto psis = esf::at_quadr<Quadr>(*psi_, is, edge);

			for (std::size_t iq = 0; iq < Quadr::size; ++iq)
			{
				// Note: eigenvectors are normalized with respect to the scalar product
				// <f | g> = <f| B |g>, for the eigenproblem A |f> = lambda B |f>.

				const auto z = (-energy + dphis[iq] + fermi_level_) / params_.temp;
				const auto f = esu::fd_int_zero(z);
				const auto fd = esu::fd_int_minus_one(z);
				const auto psi_sq = esu::sq(psis[iq]);

				density.first[iq] -= effective_dos_ * psi_sq * f;
				density.second[iq] -= effective_dos_ / params_.temp * psi_sq * fd;
			}
		}

		for (std::size_t iq = 0; iq < Quadr::size; ++iq)
			density.first[iq] += params_.dopant_conc;

		return density;
	}

	double potential_change_sup_norm() const
	{
		return esl::norm_sup(phi_.values() - prev_phi_.values());
	}

private:
	const Params& params_;
	double effective_dos_;
	double fermi_level_;

	const Poisson_solution_view phi_;
	Poisson_solution prev_phi_;

	std::optional<Schrodinger_solution_view> psi_;
};
