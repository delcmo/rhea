/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#include "JumpGradientInterface.h"

/* This function is called to compute the jump of the gradient of a given quantity when using CONTINUOUS finite element. This function acts on the sides of the cell.*/
template<>
InputParameters validParams<JumpGradientInterface>()
{
  InputParameters params = validParams<InternalSideUserObject>();

  params.addRequiredCoupledVar("variable", "the variable name this userobject is acting on.");
  params.addRequiredParam<std::string>("jump_name", "the name of the variable that will store the jump");
  
  return params;
}

JumpGradientInterface::JumpGradientInterface(const std::string & name, InputParameters parameters) :
    InternalSideUserObject(name, parameters),
    _aux(_fe_problem.getAuxiliarySystem()),
    _grad_u(coupledGradient("variable")),
    _grad_u_neighbor(coupledNeighborGradient("variable")),
    _jump_name(getParam<std::string>("jump_name")),
    _value(0.)
{
  if (_mesh.dimension()!=1)
    mooseError("The current implementation of '" << this->name() << "' can only be used with 1-D mesh.");
}

JumpGradientInterface::~JumpGradientInterface()
{
}

void
JumpGradientInterface::initialize()
{
    NumericVector<Number> & sln = _aux.solution();
    _aux.system().zero_variable(sln, _aux.getVariable(_tid, _jump_name).number());
}

void
JumpGradientInterface::execute()
{
    _value = 0.;
    NumericVector<Number> & sln = _aux.solution();

    // Compute the jump of the given variable:(grad(f_i)_x - grad(f_ip1)_x)
    for (unsigned int qp = 0; qp < _q_point.size(); ++qp) {
        _value = std::max(std::fabs(_grad_u[qp](0) - _grad_u_neighbor[qp](0)), _value);}
    
    dof_id_type _dof_nb = _current_elem->dof_number(_aux.number(), _fe_problem.getVariable(_tid, _jump_name).number(), 0);
    dof_id_type _dof_nb_neighbor = _neighbor_elem->dof_number(_aux.number(), _fe_problem.getVariable(_tid, _jump_name).number(), 0);
  
    sln.add(_dof_nb, _value);
    sln.add(_dof_nb_neighbor, _value);
}

void
JumpGradientInterface::destroy()
{
}

void
JumpGradientInterface::finalize()
{
    _aux.solution().close();
}

void
JumpGradientInterface::threadJoin(const UserObject & uo)
{
}