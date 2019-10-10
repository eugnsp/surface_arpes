#pragma once
#include <esf/boundary_cond.hpp>
#include <esf/dof/dof_mapper.hpp>
#include <esf/element/lagrange.hpp>
#include <esf/system.hpp>
#include <esf/var.hpp>
#include <esf/var_list.hpp>

#include <string>

using Schrodinger_element = esf::Lagrange<1, 1>;
using Schrodinger_dirichlet = esf::Uniform_boundary_cond<Schrodinger_element>;
using Schrodinger_var = esf::Var<Schrodinger_element, 1, Schrodinger_dirichlet, Schrodinger_dirichlet>;

class Schrodinger_system final : public esf::System<esf::Var_list<Schrodinger_var>, esf::Dof_mapper>
{
private:
	using Base = esf::System<esf::Var_list<Schrodinger_var>, esf::Dof_mapper>;

public:
	using Base::Base;

	virtual std::string name() const override
	{
		return "1D Schrodinger solver";
	}
};
