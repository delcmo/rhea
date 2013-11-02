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
This function computes the dissipative terms for all of the equations. It is dimension agnostic.
 */
template<>
InputParameters validParams<RheaArtificialVisc>()
{
  InputParameters params = validParams<Kernel>();
    // Equation and diffusion names:
    params.addParam<std::string>("equation_name", "INVALID", "Name of the equation.");
    // Boolean for dissipation in radiation equation:
    params.addParam<bool>("isRadiation", false, "boolean for dissipation in radiation equation.");
    // Material variables:
    params.addRequiredCoupledVar("density", "density");
    params.addRequiredCoupledVar("internal_energy", "internal energy");
    params.addRequiredCoupledVar("velocity", "velocity");
    params.addCoupledVar("radiation", "radiation");
  return params;
}

RheaArtificialVisc::RheaArtificialVisc(const std::string & name,
                       InputParameters parameters) :
  Kernel(name, parameters),
    // Declare equation types
    _equ_name(getParam<std::string>("equation_name")),
    _equ_type("CONTINUITY, MOMENTUM, ENERGY, RADIATION, INVALID", "INVALID"),
    // Boolean:
    _isRadiation(getParam<bool>("isRadiation")),
    // Material variables:
    _rho(coupledValue("density")),
    _grad_rho(coupledGradient("density")),
    _vel(coupledValue("velocity")),
    _grad_vel(coupledGradient("velocity")),
    _grad_rhoe(coupledGradient("internal_energy")),
    // Radiation variable:
    _grad_epsilon(isCoupled("radiation") ? coupledGradient("radiation") : _grad_zero),
    // Material property: viscosity coefficient.
    _mu(getMaterialProperty<Real>("mu")),
    _kappa(getMaterialProperty<Real>("kappa"))
{
    _equ_type = _equ_name;
}

Real RheaArtificialVisc::computeQpResidual()
{
    // Determine if cell is on boundary or not and then compute a unit vector 'l=grad(norm(vel))/norm(grad(norm(vel)))':
    Real _isOnBoundary = 1.;
    if( _current_elem->node(_i) == 0 || _current_elem->node(_i) == _mesh.nNodes()-1 )
        _isOnBoundary = 1.;

    // If statement on diffusion type:
    Real _f = _kappa[_qp]*_grad_rho[_qp](0);
    Real _g = _mu[_qp]*_rho[_qp]*_grad_vel[_qp](0);
    Real _vel2 = _vel[_qp]*_vel[_qp];
    switch (_equ_type) {
        case CONTINUITY:
            return _isOnBoundary*_f*_grad_test[_i][_qp](0);
            break;
        case MOMENTUM:
            return _isOnBoundary*(_g+_vel[_qp]*_f)*_grad_test[_i][_qp](0);
            break;
        case ENERGY:
            return _isOnBoundary*(_kappa[_qp]*_grad_rhoe[_qp](0) + 0.5*_vel2*_f + _vel[_qp]*_g)*_grad_test[_i][_qp](0);
            break;
        case RADIATION:
            return (double)_isRadiation*_isOnBoundary*_kappa[_qp]*_grad_epsilon[_qp](0)*_grad_test[_i][_qp](0);
            break;
        default:
            mooseError("INVALID equation name.");
    }
}

Real RheaArtificialVisc::computeQpJacobian()
{
  return 0.;
}

Real RheaArtificialVisc::computeQpOffDiagJacobian( unsigned int _jvar)
{
  return 0.*_jvar;
}
