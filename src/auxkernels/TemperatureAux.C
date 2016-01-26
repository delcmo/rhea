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

#include "TemperatureAux.h"
#include "MooseMesh.h"

template<>
InputParameters validParams<TemperatureAux>()
{
  InputParameters params = validParams<AuxKernel>();

  // Coupled values
  params.addRequiredCoupledVar("rho", "density");
  params.addRequiredCoupledVar("rhou", "momentum");
  params.addRequiredCoupledVar("rhoE", "energy");
  // Equation of state
  params.addRequiredParam<UserObjectName>("eos", "The parameters for the eos.");

  return params;
}

TemperatureAux::TemperatureAux(const InputParameters & parameters) :
    AuxKernel(parameters),
    // Coupled variables
    _rho(coupledValue("rho")),
    _rhou(coupledValue("rhou")),
    _rhoE(coupledValue("rhoE")),
    // User Objects for eos
    _eos(getUserObject<EquationOfState>("eos"))
{
  if (_mesh.dimension()!=1)
    mooseError("The current implementation of '" << this->name() << "' can only be used with 1-D mesh.");
}

/** Overwrite the function computevalue() */
Real
TemperatureAux::computeValue()
{
  // Compute pressure
  Real vel = _rhou[_qp]/_rho[_qp];
  Real pressure = _eos.pressure(_rho[_qp], vel, _rhoE[_qp]);

  // Return
  return _eos.temperature_from_p_rho(pressure, _rho[_qp]);
}