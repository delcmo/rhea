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

#include "RheaArtificialVisc.h"
/**
This function computes the dissipative terms for all of the equations. It only works in 1-D
 */
template<>
InputParameters validParams<RheaArtificialVisc>()
{
  InputParameters params = validParams<Kernel>();

  // Equation name:
  params.addParam<std::string>("equation_name", "invalid", "Name of the equation.");
  // Boolean for dissipation in radiation equation:
  params.addParam<bool>("isRadiation", false, "boolean for dissipation in radiation equation.");

  return params;
}

RheaArtificialVisc::RheaArtificialVisc(const std::string & name,
                       InputParameters parameters) :
  Kernel(name, parameters),
    // Declare equation types
    _equ_type("continuity x_momentum energy radiation invalid", getParam<std::string>("equation_name")),
    // Boolean:
    _isRadiation(getParam<bool>("isRadiation")),
    // Material properties:
    _kappa(getMaterialProperty<Real>("kappa")),
    _D(getMaterialProperty<Real>("diffusion"))
{
  if (_mesh.dimension()!=1)
    mooseError("The current implementation of '" << this->name() << "' can only be used with 1-D mesh.");
}

Real RheaArtificialVisc::computeQpResidual()
{
  switch (_equ_type)
  {
  case continuity:
    return _kappa[_qp]*_grad_u[_qp]*_grad_test[_i][_qp];
    break;
  case x_momentum:
    return _kappa[_qp]*_grad_u[_qp]*_grad_test[_i][_qp];
    break;
  case energy:
    return _kappa[_qp]*_grad_u[_qp]*_grad_test[_i][_qp];
    break;
  case radiation:
    if (_isRadiation && _kappa[_qp]>_D[_qp])
      return _kappa[_qp]*_grad_u[_qp]*_grad_test[_i][_qp];
    else
      return 0.;
    break;
  default:
      mooseError("INVALID equation name.");
  }
}

Real RheaArtificialVisc::computeQpJacobian()
{
  return _kappa[_qp]*_grad_phi[_j][_qp](0)*_grad_test[_i][_qp](0);
}

Real RheaArtificialVisc::computeQpOffDiagJacobian( unsigned int _jvar)
{
  return 0.*_jvar;
}