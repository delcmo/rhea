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

#include "RheaEnergy.h"
#include "MooseMesh.h"

/**
This Kernel computes the convection flux of the continuity equation.
*/
template<>
InputParameters validParams<RheaEnergy>()
{
  InputParameters params = validParams<Kernel>();

  // Boolean
  params.addParam<bool>("is_dimensional_form", true, "boolean to solve the momentum equation in a dimensional form");
  // Coupled values:
  params.addRequiredCoupledVar("rhou", "fluid momentum");
  params.addRequiredCoupledVar("rho", "fluid density");
  params.addCoupledVar("radiation", "radiation");  
  // Equation of state
  params.addRequiredParam<UserObjectName>("eos", "Equation of state");
  // Userobject computing the ICs
  params.addRequiredParam<UserObjectName>("ics", "parameters for ics.");

  return params;
}

RheaEnergy::RheaEnergy(const InputParameters & parameters) :
  Kernel(parameters),
    // Boolean
    _is_dmsl_form(getParam<bool>("is_dimensional_form")),
    // Coupled values:
    _rho(coupledValue("rho")),
    _rhou(coupledValue("rhou")),
    _epsilon(isCoupled("radiation") ? coupledValue("radiation") : _zero),
    _grad_eps(isCoupled("radiation") ?  coupledGradient("radiation") : _grad_zero),
    // Equation of state
    _eos(getUserObject<EquationOfState>("eos")),
    // Material property:
    _sigma_a(getMaterialProperty<Real>("sigma_a")),
    // Userobject computing the ICs
    _ics(getUserObject<InputFileSpecifiedICsRadHydro>("ics")),
    // Dimensional numbers
    _c(_is_dmsl_form ? _ics.c() : 1.),
    _a(_is_dmsl_form ? _ics.a() : 1.),
    // Non-dimensional numbers
    _Po(_is_dmsl_form ? 1. : _ics.P()),
    _SIGMA(_is_dmsl_form ? 1. : _ics.SIGMA_A()),
    // Integers for jacobian terms
    _rho_nb(coupled("rho")),
    _rhou_nb(coupled("rhou")),
    _epsilon_nb(coupled("radiation"))
{
  if (_mesh.dimension()!=1)
    mooseError("The current implementation of '" << this->name() << "' can only be used with 1-D mesh.");
}

Real RheaEnergy::computeQpResidual()
{
  // Convection term:
  Real vel = _rhou[_qp] / _rho[_qp];
  Real pressure = _eos.pressure(_rho[_qp], vel, _u[_qp]);
  Real conv = vel*(_u[_qp]+pressure);

  // Relaxation term:
  Real temp = _eos.temperature(_rho[_qp], vel, _u[_qp]);
  Real temp4 = temp*temp*temp*temp;
  Real relaxation = _sigma_a[_qp]*_c*(_a*temp4-_epsilon[_qp]);
  relaxation *= _SIGMA;

  // Return the flux:
  return -conv*_grad_test[_i][_qp](0) + _Po*(relaxation+vel*_grad_eps[_qp](0)/3)*_test[_i][_qp];
}

Real RheaEnergy::computeQpJacobian()
{
  // Convection term:
  Real vel = _rhou[_qp] / _rho[_qp];
  Real pressure = _eos.dp_drhoE(_rho[_qp], vel, _u[_qp]);
  Real conv = _phi[_j][_qp]*vel*(1.+pressure);

  // Relaxation term:
  Real temp = _eos.temperature(_rho[_qp], vel, _u[_qp]);
  Real temp4 = 4*temp*temp*temp*_eos.dT_drhoE(_rho[_qp], vel, _u[_qp]);
  Real relaxation = _phi[_j][_qp]*_sigma_a[_qp]*_c*_a*temp4;
  relaxation *= _SIGMA;  

  // Return jacobian term
  return -conv*_grad_test[_i][_qp](0) + relaxation*_test[_i][_qp];
}

Real RheaEnergy::computeQpOffDiagJacobian( unsigned int _jvar)
{
  if (_jvar == _rho_nb)
  {
    // Convection term:
    Real vel = _rhou[_qp] / _rho[_qp];
    Real pressure = _eos.pressure(_rho[_qp], vel, _u[_qp]);
    Real conv = -vel/_rho[_qp]*(_u[_qp]+pressure);
    pressure = _eos.dp_drho(_rho[_qp], vel, _u[_qp]);
    conv += vel*pressure;
    conv *= _phi[_j][_qp];

    // Relaxation term:
    Real temp = _eos.temperature(_rho[_qp], vel, _u[_qp]);
    Real temp4 = 4*temp*temp*temp*_eos.dT_drho(_rho[_qp], vel, _u[_qp]);
    Real relaxation = _phi[_j][_qp]*_sigma_a[_qp]*_c*_a*temp4;
    relaxation *= _SIGMA;    

    // Return jacobian term
    return -conv*_grad_test[_i][_qp](0) + (relaxation - _Po*_phi[_j][_qp]*vel/_rho[_qp]*_grad_eps[_qp](0)/3. )*_test[_i][_qp];
  }
  else if (_jvar == _rhou_nb)
  {
    // Convection term:
    Real vel = _rhou[_qp] / _rho[_qp];
    Real pressure = _eos.pressure(_rho[_qp], vel, _u[_qp]);
    Real conv = 1./_rho[_qp]*(_u[_qp]+pressure);
    pressure = _eos.dp_drhou(_rho[_qp], vel, _u[_qp]);
    conv += vel*pressure;
    conv *= _phi[_j][_qp];

    // Relaxation term:
    Real temp = _eos.temperature(_rho[_qp], vel, _u[_qp]);
    Real temp4 = 4*temp*temp*temp*_eos.dT_drhou(_rho[_qp], vel, _u[_qp]);
    Real relaxation = _phi[_j][_qp]*_sigma_a[_qp]*_c*_a*temp4;
    relaxation *= _SIGMA;

    // Return jacobian term
    return -conv*_grad_test[_i][_qp](0) + (relaxation + _Po*_phi[_j][_qp]/_rho[_qp]*_grad_eps[_qp](0)/3. )*_test[_i][_qp];
  }
  else if (_jvar == _epsilon_nb)
  {
    // Convection term:
    Real vel = _rhou[_qp] / _rho[_qp];    
    Real conv = 0.;

    // Relaxation term:
    Real relaxation = -_phi[_j][_qp]*_sigma_a[_qp]*_c*_a;
    relaxation *= _SIGMA;    

    // Return jacobian term
    return -conv*_grad_test[_i][_qp](0) + (relaxation + _Po*vel*_grad_phi[_j][_qp](0)/3. )*_test[_i][_qp];
  }
  else
    return 0.;
}