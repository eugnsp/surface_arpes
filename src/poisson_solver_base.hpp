#pragma once
#include "params.hpp"
#include "poisson_system.hpp"

#include <esf/dof/tools.hpp>
#include <esf/geometry.hpp>
#include <esf/math.hpp>
#include <esf/matrix_based/nonlinear_solver.hpp>
#include <esf/mesh/mesh1.hpp>
#include <esf/var_list.hpp>

#include <esl/dense.hpp>
#include <esl/sparse.hpp>
#include <esl/sparse/solver/pardiso_solver.hpp>
#include <esu/numeric.hpp>
#include <esu/phys.hpp>

#include <iomanip>
#include <iostream>

using Poisson_sp_solver = esl::Pardiso_solver<esl::Csr_matrix<double, esl::Symmetric_upper>>;

class Poisson_solver_base : public esf::Matrix_based_nonlinear_solver<Poisson_system, Poisson_sp_solver>
{
private:
	using Base = esf::Matrix_based_nonlinear_solver<Poisson_system, Poisson_sp_solver>;

public:
	using Base::mesh;
	using Base::system;

public:
	Poisson_solver_base(const esf::Mesh1& mesh, const Params& params) : Base(mesh), p_(params)
	{
		system().variable().set_name("phi");
	}

	void init()
	{
		system().variable().set_bnd_cond<0>(mesh(), esf::Point1{0});

		Base::init();
		esf::compute_and_set_sparsity_pattern(system(), matrix_);
	}

private:
	virtual void set_bnd_values() override
	{
		system().variable().for_each_ess_bnd_cond([this](const auto& bc) {
			for (auto& vertex : bc.vertices())
			{
				typename Base::System::template Var_vertex_dofs<0> vertex_dofs;
				system().dof_mapper().template vertex_dofs<0>(vertex, vertex_dofs);

				for (esf::Local_index i = 0; i < vertex_dofs.size(); ++i)
				{
					assert(!vertex_dofs[i].is_free);
					solution_[vertex_dofs[i].index] = bc.value();
				}
			}
		});
	}

	virtual void before_solve() override
	{
		system().variable().bnd_cond<0>().set_value(-p_.ec_surf);

		std::cout << system().name() << " non-linear solver\n"
					 "-------------------------------------\n"
					 "#    |residual|     |step|/|solution|\n"
					 "-------------------------------------"
				  << std::endl;
	}

	virtual void after_step(unsigned int n_iter, double residual_norm, double step_norm, double solution_norm) override
	{
		std::cout << std::setprecision(4) << std::setw(5) << std::left << n_iter + 1 << std::setw(15) << residual_norm
				  << step_norm / solution_norm << std::endl;
	}

protected:
	using Base::matrix_;
	using Base::solution_;

	const Params& p_;
};
