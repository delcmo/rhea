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
This function compute the Mach number. It is dimension agnostic.
**/
#include "MachNumberAux.h"

template<>
InputParameters validParams<MachNumberAux>()
{
  InputParameters params = validParams<AuxKernel>();
    // Conservative coupled variables:
    params.addRequiredCoupledVar("rho", "density");
    params.addRequiredCoupledVar("rhou", "momentum");
    // Aux coupled variables:
    params.addRequiredCoupledVar("pressure", "pressure");
    params.addRequiredParam<UserObjectName>("eos", "Equation of state");
  return params;
}

MachNumberAux::MachNumberAux(const std::string & name, InputParameters parameters) :
    AuxKernel(name, parameters),
    // Coupled variables
    _rho(coupledValue("rho")),
    _rhou(coupledValue("rhou")),
    // Aux coupled variables:
    _pressure(coupledValue("pressure")),
    // User Objects for eos
    _eos(getUserObject<EquationOfState>("eos"))
{}

Real
MachNumberAux::computeValue()
{
    // Compute the velocity:
    Real _vel = _rhou[_qp] / _rho[_qp];
    
    // Computes the speed of sounds:
    Real _c2 = _eos.c2_from_p_rho(_rho[_qp], _pressure[_qp]);

    // Return the value of the Mach number:
    return std::fabs(_vel) / std::sqrt(_c2);
}
