#pragma once
#include <es_fe/boundary_cond.hpp>
#include <es_fe/dof/dof_mapper.hpp>
#include <es_fe/element/lagrange.hpp>
#include <es_fe/system.hpp>
#include <es_fe/var.hpp>
#include <es_fe/var_list.hpp>

#include <string>

using Poisson_element = es_fe::Lagrange<1, 1>;
using Poisson_dirichlet = es_fe::Uniform_boundary_cond<Poisson_element>;
using Poisson_var = es_fe::Var<Poisson_element, 1, Poisson_dirichlet>;

class Poisson_system final : public es_fe::System<es_fe::Var_list<Poisson_var>, es_fe::Dof_mapper>
{
private:
	using Base = es_fe::System<es_fe::Var_list<Poisson_var>, es_fe::Dof_mapper>;

public:
	using Base::Base;

	virtual std::string name() const override
	{
		return "1D Poisson equation";
	}
};
