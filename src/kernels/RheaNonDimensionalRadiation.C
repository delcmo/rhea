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

#include "RheaNonDimensionalRadiation.h"
#include "MooseMesh.h"

/**
This Kernel computes the convection flux of the radiation diffusion equation
*/
template<>
InputParameters validParams<RheaNonDimensionalRadiation>()
{
  InputParameters params = validParams<Kernel>();

  // Coupled values:
  params.addRequiredCoupledVar("rho", "fluid density");  
  params.addRequiredCoupledVar("rhou", "fluid momentum");
  params.addRequiredCoupledVar("rhoE", "fluid total energy");  
  // Equation of state
  params.addRequiredParam<UserObjectName>("eos", "Equation of state");
  // Userobject computing the ICs
  params.addRequiredParam<UserObjectName>("ics", "parameters for ics.");

  return params;
}

RheaNonDimensionalRadiation::RheaNonDimensionalRadiation(const InputParameters & parameters) :
  Kernel(parameters),
    // Coupled values:
    _rho(coupledValue("rho")),
    _rhou(coupledValue("rhou")),
    _rhoE(coupledValue("rhoE")),
    // Equation of state
    _eos(getUserObject<EquationOfState>("eos")),
    // Material property:
    _sigma_a(getMaterialProperty<Real>("sigma_a")),
    _D(getMaterialProperty<Real>("diffusion")),
    // Userobject computing the ICs
    _ics(getUserObject<InputFileSpecifiedICsRadHydro>("ics")),
    // Non-dimensional numbers
    _SIGMA(_ics.SIGMA_A()),
    _K(_ics.K()),
    // Integers for jacobian terms
    _rho_nb(coupled("rho")),
    _rhou_nb(coupled("rhou")),
    _rhoE_nb(coupled("rhoE"))
{
  if (_mesh.dimension()!=1)
    mooseError("The current implementation of '" << this->name() << "' can only be used with 1-D mesh.");
}

Real RheaNonDimensionalRadiation::computeQpResidual()
{
  // Convection term:
  Real vel = _rhou[_qp]/_rho[_qp];
  Real conv = 4*vel*_u[_qp]/3.;

  // Diffusion term:
  Real diff = _D[_qp]*_grad_u[_qp](0);
  diff *= _K;

  // Relaxation term:
  Real temp = _eos.temperature(_rho[_qp], vel, _rhoE[_qp]);
  Real temp4 = temp*temp*temp*temp;
  Real relaxation = _sigma_a[_qp] * (temp4-_u[_qp]);
  relaxation *= _SIGMA;

  // Return the total expression for the radiation equation:
  return (-conv+diff)*_grad_test[_i][_qp](0) - (relaxation+vel*_grad_u[_qp](0)/3)*_test[_i][_qp];
}

Real RheaNonDimensionalRadiation::computeQpJacobian()
{
//  // Convection term:
//  Real vel = _rhou[_qp]/_rho[_qp];
//  Real conv = 4*vel*_phi[_j][_qp]/3.;
//
//  // Diffusion term:
//  Real diff = _D[_qp]*_grad_phi[_j][_qp](0);
//  diff *= _K;  
//
//  // Relaxation term:
//  Real relaxation = -_sigma_a[_qp]*_c*_phi[_j][_qp];
//  relaxation *= _SIGMA;
//
//  // Return jacobian contribution
//  return (-conv+diff)*_grad_test[_i][_qp](0) - (relaxation+vel*_grad_phi[_j][_qp](0)/3)*_test[_i][_qp];
  return 0.;
}

Real RheaNonDimensionalRadiation::computeQpOffDiagJacobian( unsigned int _jvar)
{
//  // Diffusion term:
//  Real diff = 0.;
//
//  if (_jvar == _rho_nb) // material density
//  {
//    // Convection term:
//    Real vel = _rhou[_qp]/_rho[_qp];
//    Real conv = -4*_phi[_j][_qp]*vel/_rho[_qp]*_u[_qp]/3.;
//
//    // Relaxation term:
//    Real temp = _eos.temperature(_rho[_qp], vel, _rhoE[_qp]);
//    Real temp4 = 4*_eos.dT_drho(_rho[_qp], vel, _rhoE[_qp])*temp*temp*temp;
//    Real relaxation = _sigma_a[_qp] * _c * _a*_phi[_j][_qp]*temp4;
//    relaxation *= _SIGMA;    
//
//    // Return jacobian term
//    return (-conv+diff)*_grad_test[_i][_qp](0) - (relaxation-_phi[_j][_qp]*vel/_rho[_qp]*_grad_u[_qp](0)/3)*_test[_i][_qp];
//  }
//  else if (_jvar == _rhou_nb) // material momentum
//  {
//    // Convection term:
//    Real vel = _rhou[_qp]/_rho[_qp];
//    Real conv = 4*_phi[_j][_qp]/_rho[_qp]*_u[_qp]/3.;
//
//    // Relaxation term:
//    Real temp = _eos.temperature(_rho[_qp], vel, _rhoE[_qp]);
//    Real temp4 = 4*_eos.dT_drhou(_rho[_qp], vel, _rhoE[_qp])*temp*temp*temp;
//    Real relaxation = _sigma_a[_qp] * _c * _a*_phi[_j][_qp]*temp4;
//
//    // Return jacobian term
//    return (-conv+diff)*_grad_test[_i][_qp](0) - (relaxation+_phi[_j][_qp]/_rho[_qp]*_grad_u[_qp](0)/3)*_test[_i][_qp];
//  }
//  else if (_jvar == _rhoE_nb) // material total energy
//  {
//    // Convection term:
//    Real vel = _rhou[_qp]/_rho[_qp];    
//    Real conv = 0.;
//
//    // Relaxation term:
//    Real temp = _eos.temperature(_rho[_qp], vel, _rhoE[_qp]);
//    Real temp4 = 4*_eos.dT_drhoE(_rho[_qp], vel, _rhoE[_qp])*temp*temp*temp;
//    Real relaxation = _sigma_a[_qp] * _c * _a*_phi[_j][_qp]*temp4;
//    relaxation *= _SIGMA;    
//
//    // Return jacobian term
//    return (-conv+diff)*_grad_test[_i][_qp](0) - relaxation*_test[_i][_qp];
//  }
//  else
//    return 0.;
  return 0.;
}