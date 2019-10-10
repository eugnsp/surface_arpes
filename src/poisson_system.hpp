#pragma once
#include <esf/boundary_cond.hpp>
#include <esf/dof/dof_mapper.hpp>
#include <esf/element/lagrange.hpp>
#include <esf/system.hpp>
#include <esf/var.hpp>
#include <esf/var_list.hpp>

#include <string>

using Poisson_element = esf::Lagrange<1, 1>;
using Poisson_dirichlet = esf::Uniform_boundary_cond<Poisson_element>;
using Poisson_var = esf::Var<Poisson_element, 1, Poisson_dirichlet>;

class Poisson_system final : public esf::System<esf::Var_list<Poisson_var>, esf::Dof_mapper>
{
private:
	using Base = esf::System<esf::Var_list<Poisson_var>, esf::Dof_mapper>;

public:
	using Base::Base;

	virtual std::string name() const override
	{
		return "1D Poisson equation";
	}
};
