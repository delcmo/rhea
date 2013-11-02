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

/**
This Kernel computes the convection flux of the continuity equation.
*/
template<>
InputParameters validParams<RheaEnergy>()
{
  InputParameters params = validParams<Kernel>();
    // Material values:
    params.addRequiredCoupledVar("velocity", "velocity");
    params.addCoupledVar("temperature", "temperature");
    params.addRequiredCoupledVar("pressure", "pressure");
    // Radiation value:
    params.addCoupledVar("radiation", "radiation");
    // Constant:
    params.addRequiredParam<Real>("speed_of_light", "speed of light");
    params.addRequiredParam<Real>("a", "a");
  return params;
}

RheaEnergy::RheaEnergy(const std::string & name,
                       InputParameters parameters) :
  Kernel(name, parameters),
    // Material values:
    _vel(coupledValue("velocity")),
    _temp(isCoupled("temperature") ? coupledValue("temperature") : _zero),
    _pressure(coupledValue("pressure")),
    // Radiation value:
    _epsilon(isCoupled("radiation") ? coupledValue("radiation") : _zero),
    _grad_eps(isCoupled("radiation") ?  coupledGradient("radiation") : _grad_zero),
    // Material property:
    _sigma_a(getMaterialProperty<Real>("sigma_a")),
    // Constant:
    _c(getParam<Real>("speed_of_light")),
    _a(getParam<Real>("a"))
{}

Real RheaEnergy::computeQpResidual()
{
    // Convection term:
    Real _conv = _vel[_qp]*(_u[_qp]+_pressure[_qp]);
    
    // Relaxation term:
    Real _temp4 = _temp[_qp]*_temp[_qp]*_temp[_qp]*_temp[_qp];
    Real _rel = _sigma_a[_qp] * _c * (_a*_temp4-_epsilon[_qp]);
    
    // Return the total expression for the continuity equation:
    return -_conv*_grad_test[_i][_qp](0) + (_rel+_vel[_qp]*_grad_eps[_qp](0)/3)*_test[_i][_qp];
}

Real RheaEnergy::computeQpJacobian()
{
  return ( 0 );
}

Real RheaEnergy::computeQpOffDiagJacobian( unsigned int _jvar)
{ 
    return ( 0 );
}
