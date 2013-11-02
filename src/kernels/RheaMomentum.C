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
    // Material values:
    params.addRequiredCoupledVar("velocity", "velocity");
    params.addRequiredCoupledVar("pressure", "pressure");
    // Radiation value:
    params.addCoupledVar("radiation", "radiation");
  return params;
}

RheaMomentum::RheaMomentum(const std::string & name,
                       InputParameters parameters) :
  Kernel(name, parameters),
    // Material values:
    _vel(coupledValue("velocity")),
    _pressure(coupledValue("pressure")),
    // Radiation value:
    _epsilon(isCoupled("radiation") ? coupledValue("radiation") : _zero)
{}

Real RheaMomentum::computeQpResidual()
{
    // Return the total expression for the continuity equation:
    return -( _u[_qp]*_vel[_qp] + _pressure[_qp] + _epsilon[_qp]/3 ) * _grad_test[_i][_qp](0);
}

Real RheaMomentum::computeQpJacobian()
{
  return ( 0 );
}

Real RheaMomentum::computeQpOffDiagJacobian( unsigned int _jvar)
{ 
    return ( 0 );
}
