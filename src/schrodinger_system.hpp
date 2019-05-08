#pragma once
#include <es_fe/boundary_cond.hpp>
#include <es_fe/dof/dof_mapper.hpp>
#include <es_fe/element/lagrange.hpp>
#include <es_fe/system.hpp>
#include <es_fe/var.hpp>
#include <es_fe/var_list.hpp>

#include <string>

using Schrodinger_element = es_fe::Lagrange<1, 1>;
using Schrodinger_dirichlet = es_fe::Uniform_boundary_cond<Schrodinger_element>;
using Schrodinger_var = es_fe::Var<Schrodinger_element, 1, Schrodinger_dirichlet, Schrodinger_dirichlet>;

class Schrodinger_system final : public es_fe::System<es_fe::Var_list<Schrodinger_var>, es_fe::Dof_mapper>
{
private:
	using Base = es_fe::System<es_fe::Var_list<Schrodinger_var>, es_fe::Dof_mapper>;

public:
	using Base::Base;

	virtual std::string name() const override
	{
		return "1D Schrodinger solver";
	}
};
