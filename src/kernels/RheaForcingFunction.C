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

#include "RheaForcingFunction.h"
/**
This function computes the dissipative terms for all of the equations. It is dimension agnostic.
 */
template<>
InputParameters validParams<RheaForcingFunction>()
{
  InputParameters params = validParams<Kernel>();
    // Equation and diffusion names:
    params.addParam<std::string>("equation_name", "INVALID", "Name of the equation.");
    // Constants
    params.addRequiredParam<Real>("speed_of_light", "value of the speed of light.");
    params.addRequiredParam<Real>("a", "radiation constant.");
    params.addRequiredParam<Real>("sigma_t0", "total cross-section.");
    // Equation of state:
    params.addRequiredParam<UserObjectName>("eos", "Equation of state");
  return params;
}

RheaForcingFunction::RheaForcingFunction(const std::string & name,
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
    _eos(getUserObject<EquationOfState>("eos")),
    // Material property: viscosity coefficient.
    _sigma_a(getMaterialProperty<Real>("sigma_a")),
    _D(getMaterialProperty<Real>("diffusion")),
    _mu(getMaterialProperty<Real>("mu"))
{
    _equ_type = _equ_name;
}

Real RheaForcingFunction::computeQpResidual()
{
    // Initialize some variables in order to make the code easier to read:
    Real _cos = std::cos(_q_point[_qp](0) - _t);
    Real _sin = std::sin(_q_point[_qp](0) - _t);
    //std::cout<<_q_point[_qp](0)<<std::endl;
    
    // Compute values of the pressure, temperature, ...
    Real _rho = _sin+2;
    Real _vel = _cos+2;
    Real _temp = 0.5*_eos.gamma()*(_cos+2)/(_sin+2);
    Real _temp3 = _temp*_temp*_temp;
    Real _press = _rho*_temp;
    Real _rhoe = _press / (_eos.gamma()-1);
    Real _rhoE = _rhoe + 0.5*_rho*_vel*_vel;
    Real _epsilon = _a*_temp3*_temp;
    
    // Temporal derivatives:
    Real _drhodt = -_cos;
    Real _drhoveldt = -_cos*_vel + _sin*_rho;
    Real _drhoEdt = 0.5*_eos.gamma()/(_eos.gamma()-1)*_sin + 0.5*( _vel*_drhoveldt + _rho*_vel*_sin );
    Real _depsdt = 2*_a*_temp3*_eos.gamma()*(_rho*_sin+_vel*_cos)/(_rho*_rho);
    
    // First order patial derivatives:
    Real _drhodx = _cos;
    Real _dveldx = -_sin;
    Real _dpressdx = -0.5*_eos.gamma()*_sin;
    Real _drhovel2dx = _vel*_vel*_cos - 2*_rho*_vel*_sin;
    Real _drhoEdx = -0.5*_eos.gamma()/(_eos.gamma()-1)*_sin + 0.5*_drhovel2dx;
    Real _dtempdx = 0.5*_eos.gamma()*(-_rho*_sin-_vel*_cos)/(_rho*_rho);
    Real _depsdx = 4*_a*_temp3*_dtempdx;
    
    // Second order spatial derivative:
    Real _drhodx2 = -_sin;
    Real _drhoudx2 = -(_rho*_cos+2*_sin*_cos+_vel*_sin);
    Real _drhoEdx2 = -0.5*_eos.gamma()/(_eos.gamma()-1)*_cos+0.5*(_vel*_drhoudx2 - _sin*(_vel*_drhodx + _rho*_dveldx) - _rho*_vel*_cos);
    Real _dtempdx2 = 0.5*_eos.gamma()*( (-_rho*_cos+_vel*_sin) - 4*_rho*_dtempdx*_drhodx/_eos.gamma() )/(_rho*_rho);
    Real _depsdx2 = 12*_a*_temp*_temp*_dtempdx*_dtempdx + 4*_a*_temp3*_dtempdx2;
    //std::cout<<_D[_qp]<<std::endl;
    // Source terms:
    Real _relaxation = _sigma_a[_qp]*_c*(_a*_temp3*_temp - _epsilon);
    Real _diff = _c / (3*_sigma_t);
    //std::cout<<"diff="<<_D[_qp]<<std::endl;
    //std::cout<<"depsdt="<<_depsdt<<std::endl;
    //std::cout<<"_depsdex2"<<_depsdx2<<std::endl;
    
    switch (_equ_type) {
        case CONTINUITY: /* drhodt + drhoudx */
            return -(_drhodt + _vel*_drhodx + _rho*_dveldx - _mu[_qp]*_drhodx2)*_test[_qp][_i];
            break;
        case MOMENTUM:
            return -(_drhoveldt + _drhovel2dx + _dpressdx + _depsdx/3 - _mu[_qp]*_drhoudx2)*_test[_qp][_i];
            //return -(_drhoveldt + _drhovel2dx + _dpressdx - _mu[_qp]*_drhoudx2)*_test[_qp][_i];
            break;
        case ENERGY:
            return -(_drhoEdt + (_rhoE+_press)*_dveldx + _vel*(_drhoEdx+_dpressdx) + _vel*_depsdx/3 + _relaxation- _mu[_qp]*_drhoEdx2)*_test[_qp][_i];
            //return -(_drhoEdt + (_rhoE+_press)*_dveldx + _vel*(_drhoEdx+_dpressdx) - _mu[_qp]*_drhoEdx2)*_test[_qp][_i];
            break;
        case RADIATION:
            return -(_depsdt + 4*(_vel*_depsdx+_epsilon*_dveldx)/3 - _vel*_depsdx/3 - _relaxation - _diff*_depsdx2 )*_test[_qp][_i];
            //return -(_depsdt + 4*(_vel*_depsdx+_epsilon*_dveldx)/3 - _diff*_depsdx2 )*_test[_qp][_i];
            //return -(_depsdt - _D[_qp]*_depsdx2 )*_test[_qp][_i];
            break;
        default:
            mooseError("INVALID equation name.");
    }
}

Real RheaForcingFunction::computeQpJacobian()
{
  return 0.;
}

Real RheaForcingFunction::computeQpOffDiagJacobian( unsigned int _jvar)
{
  return 0.*_jvar;
}
