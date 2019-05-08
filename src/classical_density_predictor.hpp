#pragma once
#include "params.hpp"
#include "poisson_solver_base.hpp"

#include <es_fe/mesh/mesh1.hpp>
#include <es_la/dense.hpp>
#include <es_util/numeric.hpp>

#include <cstddef>
#include <utility>

namespace es_la
{
template<typename Value, std::size_t size>
using Vector3 = Matrix<Value, size, (std::size_t)1>;
}

struct U
{
    static constexpr int i()
    {
        return 0;
    }
};

// Calculator of electron density using the Thomas-Fermi (quasi-classical) approximation
class Classical_density_predictor
{
private:
	using Potential_view = Poisson_solver_base::Solution_view<0>;

public:
	Classical_density_predictor(const Params& params) : params_(params)
	{
		effective_dos_ = effective_3d_dos(params_);
		fermi_level_ = charge_neutral_fermi_level(params_);
	}

	void set_potential_view(Potential_view phi)
	{
		phi_ = std::move(phi);
	}

	template<class Element, class Quadr, class Dofs, class R = es_la::Vector<std::pair<double, double>, Quadr::size()>>
	R get(const Dofs& dofs, const es_fe::Mesh1::Edge_view&) const
	{
		es_la::Vector<std::pair<double, double>, Quadr::size()> density;
		for (std::size_t q = 0; q < Quadr::size(); ++q)
		{
			const auto phi = phi_.template get<Element, Quadr>(q, dofs);

			const auto z = (phi + fermi_level_) / params_.lattice_temp;
			density[q].first = params_.dopant_conc - effective_dos_ * es_util::fd_int_half(z);
			density[q].second =
				-effective_dos_ / params_.lattice_temp * es_util::fd_int_minus_half(z);
		}

		return density;
	}

	template<class Element, class Quadr, class Dofs>
	es_la::Vector<std::pair<double, double>, Quadr::size()> get2(
		const Dofs& dofs, double /* edge_length */) const
	{
		es_la::Vector<std::pair<double, double>, Quadr::size()> density;
		for (std::size_t q = 0; q < Quadr::size(); ++q)
		{
			const auto phi = phi_.template get<Element, Quadr>(q, dofs);

			const auto z = (phi + fermi_level_) / params_.lattice_temp;
			density[q].first = -effective_dos_ * es_util::fd_int_half(z);
			density[q].second =
				-effective_dos_ / params_.lattice_temp * es_util::fd_int_minus_half(z);
		}

		return density;
	}

private:
	const Params& params_;
	double effective_dos_;
	double fermi_level_;

	Potential_view phi_;
};
