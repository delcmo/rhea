/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                               */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#include "RheaMomentum.h"

/**
This Kernel computes the convection flux of the continuity equation.
*/
template<>
InputParameters validParams<RheaMomentum>()
{
  InputParameters params = validParams<Kernel>();

  // Coupled values:
  params.addRequiredCoupledVar("rho", "fluid density");
  params.addRequiredCoupledVar("rhoE", "fluid total energy");
  params.addCoupledVar("radiation", "radiation");
  // Equation of state
  params.addRequiredParam<UserObjectName>("eos", "Equation of state");
  
  return params;
}

RheaMomentum::RheaMomentum(const std::string & name,
                       InputParameters parameters) :
  Kernel(name, parameters),
    // Coupled values:
    _rho(coupledValue("rho")),
    _rhoE(coupledValue("rhoE")),
    _epsilon(isCoupled("radiation") ? coupledValue("radiation") : _zero),
    // Equation of state
    _eos(getUserObject<EquationOfState>("eos"))
{
  if (_mesh.dimension()!=1)
    mooseError("The current implementation of '" << this->name() << "' can only be used with 1-D mesh.");
}

Real RheaMomentum::computeQpResidual()
{
  // Compute velocity and pressure
  Real vel = _u[_qp]/_rho[_qp];
  Real pressure = _eos.pressure(_rho[_qp], vel, _rhoE[_qp]);

  // Return
  return -(_u[_qp]*vel + pressure + _epsilon[_qp]/3)*_grad_test[_i][_qp](0);
}

Real RheaMomentum::computeQpJacobian()
{
  return 0.;
}

Real RheaMomentum::computeQpOffDiagJacobian( unsigned int _jvar)
{ 
  return 0.;
}
