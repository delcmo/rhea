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

#include "RheaForcingFunctionStream.h"
/**
This function computes the dissipative terms for all of the equations. It is dimension agnostic.
 */
template<>
InputParameters validParams<RheaForcingFunctionStream>()
{
  InputParameters params = validParams<Kernel>();
    // Equation name:
    params.addParam<std::string>("equation_name", "INVALID", "Name of the equation.");
    // Constants
    params.addRequiredParam<Real>("speed_of_light", "value of the speed of light.");
    params.addRequiredParam<Real>("a", "radiation constant.");
    params.addRequiredParam<Real>("sigma_t0", "total cross-section.");
    // Equation of state:
    params.addRequiredParam<UserObjectName>("eos", "Equation of state");
  return params;
}

RheaForcingFunctionStream::RheaForcingFunctionStream(const std::string & name,
                       InputParameters parameters) :
  Kernel(name, parameters),
    // Declare equation types
    _equ_name(getParam<std::string>("equation_name")),
    _equ_type("CONTINUITY, MOMENTUM, ENERGY, RADIATION, INVALID", "INVALID"),
    // Constant:
    _c(getParam<Real>("speed_of_light")),
    _a(getParam<Real>("a")),
    _sigma_t(getParam<Real>("sigma_t0")),
    // Equation of state:
    _eos(getUserObject<IdealGasEquationOfState>("eos")),
    // Material property: viscosity coefficient.
    _sigma_a(getMaterialProperty<Real>("sigma_a")),
    _D(getMaterialProperty<Real>("diffusion")),
    _mu(getMaterialProperty<Real>("mu"))
{
    _equ_type = _equ_name;
}

Real RheaForcingFunctionStream::computeQpResidual()
{
    // Initialize some variables in order to make the code easier to read:
    Real _cos = std::cos(_q_point[_qp](0) - _t);
    Real _cos1000 = std::cos(_q_point[_qp](0) - 1000*_t);
    Real _sin = std::sin(_q_point[_qp](0) - _t);
    Real _sin1000 = std::sin(_q_point[_qp](0) - 1000*_t);
    
    // Compute values of the pressure, temperature, ...
    Real _rho = _sin+2;
    Real _vel = 1/(_sin+2);
    Real _temp = 0.5*_eos.gamma();
    Real _temp3 = _temp*_temp*_temp;
    Real _press = _rho*_temp;
    Real _rhoe = _press / (_eos.gamma()-1);
    Real _rhoE = _rhoe + 0.5*_rho*_vel*_vel;
    Real _epsilon = _sin1000+2;
    
    // Temporal derivatives:
    Real _drhodt = -_cos;
    Real _dveldt = _cos*(_vel*_vel);
    Real _drhoveldt = 0.;
    Real _dpressdt = -0.5*_eos.gamma()*_cos;
    Real _drhoEdt = _dpressdt/(_eos.gamma()-1) + 0.5*_vel*_rho*_dveldt;
    Real _depsdt = -1000.*_cos1000;
    
    // First order patial derivatives:
    Real _drhodx = -_drhodt;
    Real _dveldx = -_dveldt;
    Real _dpressdx = -_dpressdt;
    Real _drhovel2dx = _rho*_vel*_dveldx;
    Real _drhoEdx = -_drhoEdt;
    Real _depsdx = _cos1000;
    
    // Second order spatial derivative:
    Real _drhodx2 = -_sin;
    Real _drhoudx2 = 0;
    Real _dveldx2 = (_sin*_vel+2*_dveldx*_cos)/(_vel*_vel*_vel);
    Real _drhoEdx2 = -_dpressdx/(_eos.gamma()-1) + 0.5*_rho*_vel*_dveldx2;
    Real _depsdx2 = -_sin1000;
    
    // Source terms:
    Real _relaxation = _sigma_a[_qp]*_c*(_a*_temp3*_temp - _epsilon);
    Real _diff = _c / (3*_sigma_t);
    
    switch (_equ_type) {
        case CONTINUITY: /* drhodt + drhoudx */
            return -(_drhodt + 0. - _mu[_qp]*_drhodx2)*_test[_qp][_i];
            break;
        case MOMENTUM:
            return -(_drhoveldt + _drhovel2dx + _dpressdx + _depsdx/3 - _mu[_qp]*_drhoudx2)*_test[_qp][_i];
            break;
        case ENERGY:
            return -(_drhoEdt + (_rhoE+_press)*_dveldx + _vel*(_drhoEdx+_dpressdx) + _vel*_depsdx/3 + _relaxation- _mu[_qp]*_drhoEdx2)*_test[_qp][_i];
            break;
        case RADIATION:
            return -(_depsdt + 4*(_vel*_depsdx+_epsilon*_dveldx)/3 - _vel*_depsdx/3 - _relaxation - _diff*_depsdx2 )*_test[_qp][_i];
            break;
        default:
            mooseError("INVALID equation name.");
    }
}

Real RheaForcingFunctionStream::computeQpJacobian()
{
  return 0.;
}

Real RheaForcingFunctionStream::computeQpOffDiagJacobian( unsigned int _jvar)
{
  return 0.*_jvar;
}
