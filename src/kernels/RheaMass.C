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

#include "RheaMass.h"

/**
This Kernel computes the convection flux of the continuity equation.
*/
template<>
InputParameters validParams<RheaMass>()
{
  InputParameters params = validParams<Kernel>();
  
  params.addRequiredCoupledVar("rhou", "momentum: rho*u");

  return params;
}

RheaMass::RheaMass(const std::string & name,
                       InputParameters parameters) :
  Kernel(name, parameters),
    // Coupled aux variables
    _rhou(coupledValue("rhou"))
{
  if (_mesh.dimension()!=1)
    mooseError("The current implementation of '" << this->name() << "' can only be used with 1-D mesh.");
}

Real RheaMass::computeQpResidual()
{
  // Return
  return -_rhou[_qp]*_grad_test[_i][_qp](0);
}

Real RheaMass::computeQpJacobian()
{
  return 0.;
}

Real RheaMass::computeQpOffDiagJacobian( unsigned int _jvar)
{ 
  return 0.;
}
