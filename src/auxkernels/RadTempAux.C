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
/**
This function aims at computing the Temperature at the nodes and its gradient. This auxiliary variable is coupled
to rho_bar, m_bar and E_bar defined as the product of the usual density, momentum and energy, and the cross section
A computed by the function AreaFunction.
**/
#include "RadTempAux.h"
#include "MooseMesh.h"

template<>
InputParameters validParams<RadTempAux>()
{
  InputParameters params = validParams<AuxKernel>();
  // Coupled values
  params.addRequiredCoupledVar("radiation", "radiation");
  // Constant
  params.addParam<Real>("a", 1.372e-2, "Boltzman constant");
  return params;
}

RadTempAux::RadTempAux(const InputParameters & parameters) :
    AuxKernel(parameters),
    // Coupled variables
    _eps(coupledValue("radiation")),
    // Constant
    _a(getParam<Real>("a"))
{
  if (_mesh.dimension()!=1)
    mooseError("The current implementation of '" << this->name() << "' can only be used with 1-D mesh.");
}

Real
RadTempAux::computeValue()
{
  // Return
  return std::pow(_eps[_qp]/_a, 0.25);
}