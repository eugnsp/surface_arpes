#pragma once
#include "params.hpp"
#include "poisson_solver_base.hpp"
#include "schrodinger_solver_base.hpp"

#include <es_fe/mesh/mesh1.hpp>
#include <es_la/dense.hpp>
#include <es_util/numeric.hpp>

#include <cstddef>
#include <utility>

// Calculator of electron density using the Scrodinger equation and the perturbation theory
class Quantum_density_predictor
{
private:
	using Potential_view = Poisson_solver_base::Solution_view<0>;
	using Schrodinger_view = Schrodinger_solver_base::Solution_view<0>;

public:
	Quantum_density_predictor(const Params& params) : params_(params)
	{
		effective_dos_ = effective_2d_dos(params_);
		fermi_level_ = charge_neutral_fermi_level(params_);
	}

	void xxx()
	{
		prev_phi_solution_.resize(phi_.raw().size());
		prev_phi_solution_ = phi_.raw();
	}

	void set_potential_view(Potential_view phi)
	{
		phi_ = std::move(phi);
		prev_phi_ = phi_.clone(prev_phi_solution_);
	}

	void set_schrodinger_view(Schrodinger_view psi)
	{
		psi_ = std::move(psi);
	}

	template<class Element, class Quadr, class Dofs, class R = es_la::Vector<std::pair<double, double>, Quadr::size()>>
	R get(const Dofs& dofs, const es_fe::Mesh1::Edge_view& edge) const
	{
		es_la::Vector<std::pair<double, double>, Quadr::size()> density;
		for (std::size_t q = 0; q < Quadr::size(); ++q)
		{
			const auto phi = phi_.template get<Element, Quadr>(q, dofs);
			const auto prev_phi = prev_phi_.template get<Element, Quadr>(q, dofs);

			double d = 0;
			double dd = 0;
			for (std::size_t is = 0; is < psi_.size(); ++is)
			{
				// Note: eigenvectors are normalized with respect to the scalar product
				// <f | g> = <f| B |g>, for the eigenproblem A |f> = lambda B |f>.

				const auto energy = psi_[is];
				const auto psi = psi_.template get<Element, Quadr>(is, q, edge);

				const auto z = (-energy + phi - prev_phi + fermi_level_) / params_.lattice_temp;
				const auto f = es_util::ln_one_p_exp(z);
				const auto fd = es_util::fermi(-z);

				d += -effective_dos_ * psi * psi * f;
				dd += -effective_dos_ / params_.lattice_temp * psi * psi * fd;
			}

			density[q].first = params_.dopant_conc + d;
			density[q].second = dd;
		}

		return density;
	}

	template<class Element, class Quadr, class Dofs>
	es_la::Vector<std::pair<double, double>, Quadr::size()> get2(const Dofs& dofs, double edge_length) const
	{
		es_la::Vector<std::pair<double, double>, Quadr::size()> density;
		for (std::size_t q = 0; q < Quadr::size(); ++q)
		{
			const auto phi = phi_.template get<Element, Quadr>(q, dofs);
			const auto prev_phi = prev_phi_.template get<Element, Quadr>(q, dofs);

			double d = 0;
			double dd = 0;
			for (std::size_t is = 0; is < psi_.size(); ++is)
			{
				const auto energy = psi_[is];
				const auto psi = psi_.template get<Element, Quadr>(is, q, dofs);

				const auto z = (-energy + phi - prev_phi + fermi_level_) / params_.lattice_temp;
				const auto f = es_util::ln_one_p_exp(z);
				const auto fd = es_util::fermi(-z);

				d += -effective_dos_ * psi * psi * f;
				dd += -effective_dos_ / params_.lattice_temp * psi * psi * fd;

				//std::cout << z << '\n';
			}

			density[q].first = d;
			density[q].second = dd;
		}

		return density;
	}

private:
	const Params& params_;
	double effective_dos_;
	double fermi_level_;

	Potential_view phi_;
	Potential_view prev_phi_;
	es_la::Vector_xd prev_phi_solution_;

	Schrodinger_view psi_;
};
