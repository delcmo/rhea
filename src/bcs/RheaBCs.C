#include "RheaBCs.h"

template<>
InputParameters validParams<RheaBCs>()
{
  InputParameters params = validParams<IntegratedBC>();

    params.addRequiredParam<std::string>("equation_name", "The name of the equation this BC is acting on");
    // Coupled variables:
    params.addRequiredCoupledVar("velocity", "velocity");
    params.addRequiredCoupledVar("density", "density");
    params.addRequiredCoupledVar("pressure", "pressure");
    params.addCoupledVar("epsilon", "epsilon");
    // Input parameters:
    params.addRequiredParam<Real>("p_bc", "Static pressure at the boundary");
    params.addParam<Real>("T_bc", 0.0, "Static temperature at the boundary");
    params.addParam<Real>("v_bc", 0.0, "Velocity at the boundary.");
    params.addParam<Real>("eps_left", 0.0, "left value of the radiation");
    params.addParam<Real>("eps_right", 0.0, "right value of the radiation");
    params.addRequiredParam<Real>("speed_of_light", "speed of light");
    // Equation of state:
    params.addRequiredParam<UserObjectName>("eos", "The name of equation of state object to use.");
  return params;
}

RheaBCs::RheaBCs(const std::string & name, InputParameters parameters) :
    IntegratedBC(name, parameters),
    // Name of the equation:
    _eqn_name(getParam<std::string>("equation_name")),
    _eqn_type("CONTINUITY, MOMENTUM, ENERGY, RADIATION, INVALID", "INVALID"),
    // Coupled variables:
    _vel(coupledValue("velocity")),
    _rho(coupledValue("density")),
    _pressure(coupledValue("pressure")),
    _epsilon(isCoupled("epsilon") ? coupledValue("epsilon") : _zero),
    _vel_old(coupledValueOld("velocity")),
    _rho_old(coupledValueOld("density")),
    _pressure_old(coupledValueOld("pressure")),
    // Boundary condition parameters:
    _p_bc(getParam<Real>("p_bc")),
    _T_bc(getParam<Real>("T_bc")),
    _v_bc(getParam<Real>("v_bc")),
    _eps_left(getParam<Real>("eps_left")),
    _eps_right(getParam<Real>("eps_right")),
    _c_light(getParam<Real>("speed_of_light")),
    // Material:
    _D(getMaterialProperty<Real>("diffusion")),
    // Equation of state:
    _eos(getUserObject<EquationOfState>("eos"))
{
    _eqn_type = _eqn_name;
}

Real
RheaBCs::computeQpResidual()
{
    // Compute the mach number:
    Real _Mach = _vel[_qp] / std::sqrt(_eos.c2_from_p_rho(_rho[_qp], _pressure[_qp]));
    Real _Mach_old = _vel_old[_qp] / std::sqrt(_eos.c2_from_p_rho(_rho_old[_qp], _pressure_old[_qp]));
    //std::cout<<"Mach="<<_Mach<<std::endl;
    
    // Declare variables  that are used in the fluxes:
    Real _rho_bc = 0.; Real _rhou_bc = 0.; Real _vel_bc = 0.; Real _e_bc = 0.; Real _press_bc = 0.; Real _rhoE_bc = 0.;
    
    // inlet or outlet:
    if ( _vel[_qp]*_normals[_qp](0) <0 ) // Inlet
    {
        if (_Mach_old <= 1) { // subsonic
            _press_bc = _p_bc;
            _rho_bc = _eos.rho_from_p_T(_press_bc, _T_bc);
            _e_bc = _eos.e_from_p_rho(_press_bc, _rho_bc);
            _rhou_bc = _u[_qp];
            _rhoE_bc = _u[_qp];
            _vel_bc = _vel[_qp];
        }
        else { // supersonic
            _press_bc = _p_bc;
            _rho_bc = _eos.rho_from_p_T(_press_bc, _T_bc);
            _e_bc = _eos.e_from_p_rho(_press_bc, _rho_bc);
            _vel_bc = _v_bc;
            _rhou_bc = _rho_bc*_vel_bc;
            _rhoE_bc = _rho_bc*(_e_bc+0.5*_vel_bc*_vel_bc);
        }
    }
    else // outlet
    {
        if (_Mach_old < 1) { // subsonic
            _press_bc = _p_bc;
            _rho_bc = _u[_qp];//_rho[_qp];
            _rhou_bc = _u[_qp];//_rhou[_qp];
            _vel_bc = _vel[_qp];//_rhou[_qp]/_rho[_qp];
            _rhoE_bc = _u[_qp];//_rhoE[_qp];
        }
        else { // supersonic
            //std::cout<<"supersonic"<<std::endl;
            _press_bc = _pressure[_qp];
            _rho_bc = _u[_qp];//_rho[_qp];
            _rhou_bc = _u[_qp];//_rhou[_qp];
            _vel_bc = _vel[_qp];//_rhou[_qp]/_rho[_qp];
            _rhoE_bc = _u[_qp];//_rhoE[_qp];
        }
    }
    
    // Switch statement for type of equation:
    Real _eps = 0.;
    switch (_eqn_type) {
        case CONTINUITY:
            return _rho_bc*_vel_bc*_normals[_qp](0)*_test[_i][_qp];
            break;
        case MOMENTUM:
            return (_rhou_bc*_vel_bc +_press_bc + _epsilon[_qp]/3)*_normals[_qp](0)*_test[_i][_qp];
            break;
        case ENERGY:
            return (_rhoE_bc+_press_bc)*_vel_bc*_normals[_qp](0)*_test[_i][_qp];
            break;
        case RADIATION:
            _eps = -0.5*(_normals[_qp](0)-1)*_eps_left + 0.5*(1+_normals[_qp](0))*_eps_right;
            return (4*_vel[_qp]*_u[_qp]/3*_normals[_qp](0)+0.5*_D[_qp]*(_u[_qp]-_eps))*_test[_i][_qp];
            break;
        default:
            mooseError("The equation with name: \"" << _eqn_name << "\" is not supported in the \"RheaBCs\" type of boundary condition.");
            return 0.;
            break;
    }
}

Real
RheaBCs::computeQpJacobian()
{
  // TODO
  return 0;
}

Real
RheaBCs::computeQpOffDiagJacobian(unsigned jvar)
{
  // TODO
  return 0;
}
