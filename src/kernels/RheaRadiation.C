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
This Kernel computes the convection flux of the continuity equation.
*/
template<>
InputParameters validParams<RheaRadiation>()
{
  InputParameters params = validParams<Kernel>();
    // Material values:
    params.addRequiredCoupledVar("velocity", "velocity");
    params.addRequiredCoupledVar("temperature", "temperature");
    // Constant:
    params.addRequiredParam<Real>("speed_of_light", "speed of light");
    params.addRequiredParam<Real>("a", "a");
  return params;
}

RheaRadiation::RheaRadiation(const std::string & name,
                       InputParameters parameters) :
  Kernel(name, parameters),
    // Material values:
    _vel(coupledValue("velocity")),
    _temp(coupledValue("temperature")),
    // Material property:
    _sigma_a(getMaterialProperty<Real>("sigma_a")),
    _D(getMaterialProperty<Real>("diffusion")),
    // Constant:
    _c(getParam<Real>("speed_of_light")),
    _a(getParam<Real>("a"))
{}

Real RheaRadiation::computeQpResidual()
{
    // Convection term:
    Real _conv = 4*_vel[_qp]*_u[_qp]/3;
    
    // Diffusion term:
    Real _diff = _D[_qp]*_grad_u[_qp](0);
    //std::cout<<"diff kernel="<<_D[_qp]<<std::endl;
    
    // Relaxation term:
    Real _temp4 = _temp[_qp]*_temp[_qp]*_temp[_qp]*_temp[_qp];
    Real _rel = _sigma_a[_qp] * _c * (_a*_temp4-_u[_qp]);
    
    // Return the total expression for the continuity equation:
    return -(_conv-_diff)*_grad_test[_i][_qp](0) - (_rel+_vel[_qp]*_grad_u[_qp](0)/3)*_test[_i][_qp];
}

Real RheaRadiation::computeQpJacobian()
{
  return ( 0 );
}

Real RheaRadiation::computeQpOffDiagJacobian( unsigned int _jvar)
{ 
    return ( 0 );
}
