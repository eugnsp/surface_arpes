#pragma once
#include "params.hpp"
#include "poisson_solver_base.hpp"
#include "schrodinger_system.hpp"

#include <es_fe/dof/tools.hpp>
#include <es_fe/geometry.hpp>
#include <es_fe/io/matlab_writer1.hpp>
#include <es_fe/math.hpp>
#include <es_fe/matrix_based/eigen_solver.hpp>
#include <es_fe/mesh/mesh1.hpp>
#include <es_fe/var_list.hpp>

#include <es_la/dense.hpp>
#include <es_la/sparse.hpp>
#include <es_la/sparse/solver/feast_interval_solver.hpp>
#include <es_la/sparse/solver/feast_solver2.hpp>
#include <es_util/numeric.hpp>
#include <es_util/phys.hpp>

using Schrodinger_eigen_solver = es_la::Feast_interval_solver<es_la::Csr_matrix<double, es_la::Symmetric_upper>>;

class Schrodinger_solver_base : public es_fe::Matrix_based_eigen_solver<Schrodinger_system, Schrodinger_eigen_solver>
{
private:
	using Base = es_fe::Matrix_based_eigen_solver<Schrodinger_system, Schrodinger_eigen_solver>;

public:
	using Base::mesh;
	using Base::system;

public:
	Schrodinger_solver_base(const es_fe::Mesh1& mesh, const Params& params) : Base(mesh), p_(params)
	{
		system().variable().set_name("psi");
	}

	void init()
	{
		system().variable().set_bnd_cond<0>(mesh(), es_fe::Point1{0}, 0);
		system().variable().set_bnd_cond<1>(mesh(), es_fe::Point1{p_.length}, 0);

		Base::init();
		es_fe::compute_and_set_sparsity_pattern(system(), matrix_a_);
		es_fe::compute_and_set_sparsity_pattern(system(), matrix_b_);
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
