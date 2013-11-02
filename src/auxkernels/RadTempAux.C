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

template<>
InputParameters validParams<RadTempAux>()
{
  InputParameters params = validParams<AuxKernel>();
    params.addRequiredCoupledVar("temperature", "temperature");
    params.addRequiredParam<Real>("a", "a");
  return params;
}

RadTempAux::RadTempAux(const std::string & name, InputParameters parameters) :
    AuxKernel(name, parameters),
  // Coupled variables
   _temperature(coupledValue("temperature")),
   _a(getParam<Real>("a"))
{}

Real
RadTempAux::computeValue()
{
    return _a*std::pow(_temperature[_qp], 4);
}
