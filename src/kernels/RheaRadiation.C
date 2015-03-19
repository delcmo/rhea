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

#include "RheaRadiation.h"

/**
This Kernel computes the convection flux of the radiation diffusion equation
*/
template<>
InputParameters validParams<RheaRadiation>()
{
  InputParameters params = validParams<Kernel>();

  // Coupled values:
  params.addRequiredCoupledVar("rhou", "fluid momentum");
  params.addRequiredCoupledVar("rho", "fluid density");
  params.addRequiredCoupledVar("rhoE", "fluid total energy");  
  // Equation of state
  params.addRequiredParam<UserObjectName>("eos", "Equation of state");
  // Constant:
  params.addParam<Real>("speed_of_light", 299.792, "speed of light");
  params.addParam<Real>("a", 1.372e-2, "Boltzman constant");

  return params;
}

RheaRadiation::RheaRadiation(const std::string & name,
                       InputParameters parameters) :
  Kernel(name, parameters),
    // Coupled values:
    _rho(coupledValue("rho")),
    _rhou(coupledValue("rhou")),
    _rhoE(coupledValue("rhoE")),
    // Equation of state
    _eos(getUserObject<EquationOfState>("eos")),
    // Material property:
    _sigma_a(getMaterialProperty<Real>("sigma_a")),
    _D(getMaterialProperty<Real>("diffusion")),
    // Constant:
    _c(getParam<Real>("speed_of_light")),
    _a(getParam<Real>("a"))
{
  if (_mesh.dimension()!=1)
    mooseError("The current implementation of '" << this->name() << "' can only be used with 1-D mesh.");
}

Real RheaRadiation::computeQpResidual()
{
  // Convection term:
  Real vel = _rhou[_qp]/_rho[_qp];
  Real conv = 4*vel*_u[_qp]/3.;

  // Diffusion term:
  Real diff = _D[_qp]*_grad_u[_qp](0);

  // Relaxation term:
  Real temp = _eos.temperature(_rho[_qp], vel, _rhoE[_qp]);
  Real temp4 = temp*temp*temp*temp;
  Real relaxation = _sigma_a[_qp] * _c * (_a*temp4-_u[_qp]);

  // Return the total expression for the continuity equation:
  return (-conv+diff)*_grad_test[_i][_qp](0) - (relaxation+vel*_grad_u[_qp](0)/3)*_test[_i][_qp];
}

Real RheaRadiation::computeQpJacobian()
{
  return 0.;
}

Real RheaRadiation::computeQpOffDiagJacobian( unsigned int _jvar)
{ 
  return 0.;
}
