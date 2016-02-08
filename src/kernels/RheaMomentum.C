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
#include "MooseMesh.h"

/**
This Kernel computes the convection flux of the continuity equation.
*/
template<>
InputParameters validParams<RheaMomentum>()
{
  InputParameters params = validParams<Kernel>();

  // Boolean
  params.addParam<bool>("is_dimensional_form", true, "boolean to solve the momentum equation in a dimensional form");
  // Coupled values:
  params.addRequiredCoupledVar("rho", "fluid density");
  params.addRequiredCoupledVar("rhoE", "fluid total energy");
  params.addCoupledVar("radiation", "radiation");
  // Equation of state
  params.addRequiredParam<UserObjectName>("eos", "Equation of state");
  // Userobject computing the ICs
  params.addRequiredParam<UserObjectName>("ics", "parameters for ics.");
  
  return params;
}

RheaMomentum::RheaMomentum(const InputParameters & parameters) :
  Kernel(parameters),
    // Boolean
    _is_dmsl_form(getParam<bool>("is_dimensional_form")),
    // Coupled values:
    _rho(coupledValue("rho")),
    _rhoE(coupledValue("rhoE")),
    _epsilon(isCoupled("radiation") ? coupledValue("radiation") : _zero),
    // Equation of state
    _eos(getUserObject<EquationOfState>("eos")),
    // Userobject computing the ICs
    _ics(getUserObject<InputFileSpecifiedICsRadHydro>("ics")),
    // Non-dimensional number Po
    _Po(_is_dmsl_form ? 1. : _ics.P()),
    // Integers for jacobian terms
    _rho_nb(coupled("rho")),
    _rhoE_nb(coupled("rhoE")),
    _epsilon_nb(coupled("radiation"))
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
  return -(_u[_qp]*vel + pressure + _Po*_epsilon[_qp]/3)*_grad_test[_i][_qp](0);
}

Real RheaMomentum::computeQpJacobian()
{
  // Compute velocity and pressure
  Real vel = _u[_qp]/_rho[_qp];
  Real pressure = _eos.dp_drhou(_rho[_qp], vel, _rhoE[_qp]);

  // Return jacobian terms
  return -(2.*vel+pressure)*_phi[_j][_qp]*_grad_test[_i][_qp](0);
}

Real RheaMomentum::computeQpOffDiagJacobian( unsigned int _jvar)
{
  if (_jvar == _rho_nb) // material density
  {
    // Compute velocity and pressure
    Real vel = _u[_qp]/_rho[_qp];
    Real pressure = _eos.dp_drho(_rho[_qp], vel, _rhoE[_qp]);

    // Return jacobian terms
    return -(-vel*vel+pressure)*_phi[_j][_qp]*_grad_test[_i][_qp](0);
  }
  else if (_jvar == _rhoE_nb) // material total energy
  {
    // Compute velocity and pressure
    Real vel = _u[_qp]/_rho[_qp];
    Real pressure = _eos.dp_drhoE(_rho[_qp], vel, _rhoE[_qp]);

    // Return jacobian terms
    return -pressure*_phi[_j][_qp]*_grad_test[_i][_qp](0);
  }
  else if (_jvar == _epsilon_nb) // radiation
  {
    // Return jacobian term
    return -_phi[_j][_qp]/3.*_grad_test[_i][_qp](0);
  }
  else
    return 0.;
}