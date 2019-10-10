#pragma once
#include "params.hpp"
#include "poisson_solver_base.hpp"
#include "schrodinger_system.hpp"

#include <esf/dof/tools.hpp>
#include <esf/geometry.hpp>
#include <esf/math.hpp>
#include <esf/matrix_based/eigen_solver.hpp>
#include <esf/mesh/mesh1.hpp>

#include <esl/dense.hpp>
#include <esl/sparse.hpp>
#include <esl/sparse/solver/feast_interval_solver.hpp>

using Schrodinger_eigen_solver = esl::Feast_interval_solver<esl::Csr_matrix<double, esl::Symmetric_upper>>;

class Schrodinger_solver_base : public esf::Matrix_based_eigen_solver<Schrodinger_system, Schrodinger_eigen_solver>
{
private:
	using Base = esf::Matrix_based_eigen_solver<Schrodinger_system, Schrodinger_eigen_solver>;

public:
	using Base::mesh;
	using Base::system;

public:
	Schrodinger_solver_base(const esf::Mesh1& mesh, const Params& params) : Base(mesh), p_(params)
	{
		system().variable().set_name("psi");
	}

	void init()
	{
		system().variable().set_bnd_cond<0>(mesh(), esf::Point1{0}, 0);
		system().variable().set_bnd_cond<1>(mesh(), esf::Point1{p_.length}, 0);

		Base::init();
		esf::compute_and_set_sparsity_pattern(system(), matrix_a_);
		esf::compute_and_set_sparsity_pattern(system(), matrix_b_);
	}

	void solve()
	{
		Base::solve();
	}

protected:
	using Base::matrix_a_;
	using Base::matrix_b_;

	const Params& p_;
};
