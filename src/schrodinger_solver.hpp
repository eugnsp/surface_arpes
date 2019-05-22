#pragma once
#include "params.hpp"
#include "poisson_solver_base.hpp"
#include "schrodinger_solver_base.hpp"

#include <es_fe/dof/tools.hpp>
#include <es_fe/geometry.hpp>
#include <es_fe/io/matlab_writer1.hpp>
#include <es_fe/math.hpp>
#include <es_fe/matrix_based/eigen_solver.hpp>
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

class Schrodinger_solver final : public Schrodinger_solver_base
{
private:
	using Potential_view = Poisson_solver_base::Solution_view;

public:
	Schrodinger_solver(const es_fe::Mesh1& mesh, const Params& params, Potential_view phi) :
		Schrodinger_solver_base(mesh, params), phi_(phi)
	{}

private:
	virtual std::pair<double, double> eigen_values_range() const override
	{
		auto max_v = phi_[es_fe::Vertex_index{0}];
		for (auto& vertex : mesh().vertices())
			max_v = std::max(max_v, phi_[*vertex]);

		const auto max_e = charge_neutral_fermi_level(p_) + 7 * p_.lattice_temp;
		return {-max_v, max_e};
	}

	virtual unsigned int eigen_space_dim() const override
	{
		const auto max_e = charge_neutral_fermi_level(p_) + 7 * p_.lattice_temp;

		const auto n_states = es_util::trapez_int(
								  mesh().n_vertices(), [this](auto i) { return mesh().vertex(i).x(); },
								  [this, max_e](auto i) {
									  const auto u = max_e + phi_[i];
									  return u <= 0 ? 0. : std::sqrt(2 * p_.effective_mass * u);
								  },
								  0.) /
							  es_util::math::pi;

		return static_cast<unsigned int>(5 + 1.25 * n_states);
	}

	virtual void assemble() override
	{
		matrix_a_.zero();
		matrix_b_.zero();

		for (const auto& face : mesh().edges())
			assemble_on_edge(face);
	}

	void assemble_on_edge(const es_fe::Mesh1::Edge_view& edge)
	{
		const auto length = es_fe::length(edge);

		using Stiff_quadr = es_fe::Quadr<2 * (Schrodinger_element::order - 1), 1>;
		using Mass_quadr = es_fe::Quadr<2 * Schrodinger_element::order, 1>;

		const auto grads = es_fe::gradients<Schrodinger_element, Stiff_quadr>(es_fe::inv_jacobian(edge));
		const auto stiffness_matrix =
			es_fe::stiffness_matrix<Schrodinger_element, Stiff_quadr>(grads, length / (2 * p_.effective_mass));

		const auto dofs = system().dofs(edge);

		const auto mass_matrix = es_fe::mass_matrix<Schrodinger_element, Mass_quadr>(length);

		const auto potential_matrix = es_fe::mass_matrix<Schrodinger_element, Mass_quadr>(
			[this, &edge](auto q) { return phi_.template get<Schrodinger_element, Mass_quadr>(q, edge); }, -length);

		for (es_fe::Local_index c = 0; c < dofs.size(); ++c)
			if (dofs[c].is_free)
			{
				for (es_fe::Local_index r = 0; r <= c; ++r)
					if (dofs[r].is_free)
					{
						auto ir = dofs[r].index;
						auto ic = dofs[c].index;
						es_util::sort2(ir, ic);

						matrix_a_(ir, ic) += stiffness_matrix(r, c) + potential_matrix(r, c);
						matrix_b_(ir, ic) += mass_matrix(r, c);
					}
			}
	}

public:
	void write(const std::string& file_name)
	{
		es_la::Matfile_writer m(file_name);
		m.write("psi", eigen_vectors_);
		m.write("en", eigen_values_);
	}

private:
	const Potential_view phi_;
};
