#pragma once
#include "params.hpp"
#include "poisson_solver_base.hpp"

#include <es_fe/geometry.hpp>
#include <es_fe/io/matlab_writer1.hpp>
#include <es_fe/math.hpp>
#include <es_fe/mesh/mesh1.hpp>
#include <es_fe/var_list.hpp>

#include <es_la/dense.hpp>
#include <es_la/sparse.hpp>
#include <es_util/numeric.hpp>
#include <es_util/phys.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#include <iostream>

template<class Density_predictor>
class Poisson_solver final : public Poisson_solver_base
{
public:
	Poisson_solver(const es_fe::Mesh1& mesh, const Params& params) :
		Poisson_solver_base(mesh, params), density_predictor_(params, solution_view())
	{}

	void init()
	{
		Poisson_solver_base::init();
		solution_ = 0; // Local charge neutrality at the beginning
	}

	void solve()
	{
		density_predictor_.before_solve();
		Poisson_solver_base::solve();
	}

	Density_predictor& density_predictor()
	{
		return density_predictor_;
	}

private:
	virtual void assemble() override
	{
		matrix_.zero();
		rhs_ = 0;

		for (const auto& face : mesh().edges())
			assemble_on_edge(face);
	}

	void assemble_on_edge(const es_fe::Mesh1::Edge_view& edge)
	{
		const auto eps = this->p_.eps;
		const auto length = es_fe::length(edge);

		using Stiff_quadr = es_fe::Quadr<2 * (Poisson_element::order - 1), 1>;
		using Mass_quadr = es_fe::Quadr<2 * Poisson_element::order, 1>;

		const auto grads = es_fe::gradients<Poisson_element, Stiff_quadr>(es_fe::inv_jacobian(edge));
		const auto stiffness_matrix = es_fe::stiffness_matrix<Poisson_element, Stiff_quadr>(grads, length * eps);

		const auto dofs = system().dofs(edge);
		const auto density = density_predictor_.template get<Poisson_element, Mass_quadr>(dofs, edge);

		const auto mass_matrix = es_fe::mass_matrix<Poisson_element, Mass_quadr>(
			[&density](auto q) { return density[q].second; }, es_util::math::four_pi * length);
		const auto load_vector = es_fe::load_vector<Poisson_element, Mass_quadr>(
			[&density](auto q) { return density[q].first; }, es_util::math::four_pi * length);

		for (es_fe::Local_index c = 0; c < dofs.size(); ++c)
		{
			if (dofs[c].is_free)
			{
				for (es_fe::Local_index r = 0; r <= c; ++r)
					if (dofs[r].is_free)
					{
						auto ir = dofs[r].index;
						auto ic = dofs[c].index;
						es_util::sort2(ir, ic);

						matrix_(ir, ic) += stiffness_matrix(r, c) - mass_matrix(r, c);
					}

				rhs_[dofs[c].index] += load_vector[c];
			}

			for (es_fe::Local_index r = 0; r < dofs.size(); ++r)
				if (dofs[r].is_free)
				{
					auto ir = dofs[r].index;
					auto ic = dofs[c].index;
					rhs_[ir] += -stiffness_matrix(r, c) * solution_[ic];
				}
		}
	}

public:
	void write(const std::string& file_name)
	{
		using namespace es_util::au::literals;

		es_la::Vector_xd ec(*mesh().n_vertices(), 0);
		es_la::Vector_xd n(*mesh().n_vertices(), 0);
		// es_la::Vector_xd n(*mesh().n_edges(), 0);

		for (es_fe::Vertex_index vertex{0}; vertex < mesh().n_vertices(); ++vertex)
		{
			System::template Var_vertex_dofs<0> vertex_dofs;
			system().dof_mapper().template vertex_dofs<0>(vertex, vertex_dofs);
			ec[*vertex] = -es_util::au::to_volt(solution_[vertex_dofs[0].index]);
			// ec[*vertex] = -solution_[vertex_dofs[0].index];
			// n[*vertex] = es_util::au::to_per_cm3(density_predictor_.get2(vertex));
		}

		for (auto& edge : mesh().edges())
		{
			using Quadr = es_fe::Quadr<1, 1>;

			const auto dofs = system().dofs(edge);
			const auto density = density_predictor_.template get<Poisson_element, Quadr>(dofs, edge);

			n[**edge] = es_util::au::to_per_cm3(density[0].first);
		}

		es_fe::Matlab_writer1 m(file_name, mesh(), 1_nm);
		m.write_vertex_field("ec", ec);
		m.write_vertex_field("n", n);

		const auto fermi_level = charge_neutral_fermi_level(p_);
		//		m.write_scalar("f", es_util::au::to_volt(fermi_level));
		m.write_scalar("f", fermi_level);
	}

private:
	Density_predictor density_predictor_;
};
