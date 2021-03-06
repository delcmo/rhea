#include "RheaBCs.h"

template<>
InputParameters validParams<RheaBCs>()
{
  InputParameters params = validParams<IntegratedBC>();

  // Boolean
  params.addParam<bool>("is_dimensional_form", true, "boolean to solve the momentum equation in a dimensional form");
  // Equation name
  params.addRequiredParam<std::string>("equation_name", "The name of the equation this BC is acting on");
  // Coupled variables:
  params.addRequiredCoupledVar("rho", "density");
  params.addRequiredCoupledVar("rhou", "momentum");
  params.addRequiredCoupledVar("rhoE", "energy");  
  params.addCoupledVar("epsilon", "epsilon");
  // Coupled aux variables
  params.addRequiredCoupledVar("pressure", "pressure");
  // Equation of state:
  params.addRequiredParam<UserObjectName>("eos", "The name of equation of state object to use.");
  // Userobject computing the ICs
  params.addRequiredParam<UserObjectName>("ics", "parameters for ics.");

  return params;
}

RheaBCs::RheaBCs(const InputParameters & parameters) :
    IntegratedBC(parameters),
    // Name of the equation:
    _eqn_type("continuity x_momentum energy radiation invalid", getParam<std::string>("equation_name")),
    // Coupled variables:
    _rho(coupledValue("rho")),
    _rho_old(coupledValueOld("rho")),
    _rhou(coupledValue("rhou")),
    _rhou_old(coupledValueOld("rhou")),
    _epsilon(isCoupled("epsilon") ? coupledValue("epsilon") : _zero),
    _epsilon_old(isCoupled("epsilon") ? coupledValueOld("epsilon") : _zero),
    // Coupled aux variables
    _pressure(coupledValue("pressure")),
    _pressure_old(coupledValueOld("pressure")),
    // Material:
    _D(getMaterialProperty<Real>("diffusion")),
    // Equation of state:
    _eos(getUserObject<IdealGasEquationOfState>("eos")),
    // Userobject computing the ICs
    _ics(getUserObject<InputFileSpecifiedICsRadHydro>("ics")),
    // Non-dimensional numbers
    _Po(getParam<bool>("is_dimensional_form") ? 1. : _ics.P()),
    _K(getParam<bool>("is_dimensional_form") ? 1. : _ics.K()),
    // Integers for jacobian terms
    _rho_nb(coupled("rho")),
    _rhou_nb(coupled("rhou")),
    _rhoE_nb(coupled("rhoE")),
    _epsilon_nb(coupled("epsilon"))
{
  // Compute press_hat_pre and press_hat_post
  _press_hat_pre = _eos.p_from_T_rho(_ics.T_hat_pre(), _ics.rho_hat_pre());
  _press_hat_post = _eos.p_from_T_rho(_ics.T_hat_post(), _ics.rho_hat_post());
}

Real
RheaBCs::computeQpResidual()
{
  // Compute the mach number:
  Real vel = _rhou[_qp]/_rho[_qp];
  Real vel_old = _rhou_old[_qp]/_rho_old[_qp];
  Real Mach_old = vel_old / std::sqrt(_eos.c2_from_p_rho(_rho_old[_qp], _pressure_old[_qp], _epsilon_old[_qp]));

  // Declare variables  that are used in the fluxes:
  Real rho_bc, rhou_bc, vel_bc, e_bc, press_bc, rhoE_bc, epsilon_bc;

  // inlet or outlet:
  if ( vel*_normals[_qp](0) <0 ) // Inlet
  {
    if (Mach_old <= 1) // subsonic inlet
    {
      press_bc = _press_hat_pre;
      rho_bc = _eos.rho_from_p_T(press_bc, _ics.T_hat_pre());
      e_bc = _eos.e_from_p_rho(press_bc, rho_bc);
      rhou_bc = _rhou[_qp];
      rhoE_bc = _u[_qp];
      vel_bc = vel;
      epsilon_bc = _ics.eps_hat_pre();
    }
    else // supersonic inlet
    {
      press_bc = _press_hat_pre;
      rho_bc = _eos.rho_from_p_T(press_bc, _ics.T_hat_pre());
      e_bc = _eos.e_from_p_rho(press_bc, rho_bc);
      vel_bc = _ics.vel_hat_pre();
      rhou_bc = rho_bc*vel_bc;
      rhoE_bc = rho_bc*(e_bc+0.5*vel_bc*vel_bc);
      epsilon_bc = _ics.eps_hat_pre();
    }
  }
  else // outlet
  {
//    if (Mach_old < 1) // subsonic outlet
//    {
//      press_bc = _press_hat_post;
//      rho_bc = _rho[_qp];
//      rhou_bc = _rhou[_qp];
//      vel_bc = vel;
//      rhoE_bc = _u[_qp];
//      epsilon_bc = _epsilon[_qp];
//    }
//    else // supersonic outlet
//    {

//      press_bc = _press_hat_post; // _pressure[_qp];
//      rho_bc = _rho[_qp];
//      rhou_bc = _rhou[_qp];
//      vel_bc = vel;
//      rhoE_bc = _u[_qp];
//      epsilon_bc = _epsilon[_qp];
    press_bc = _pressure[_qp];
    rho_bc = _rho[_qp];
    vel_bc = _ics.mach_hat_post()*std::sqrt(_eos.c2_from_p_rho(rho_bc, press_bc, 0.));
    rhou_bc = rho_bc*vel_bc;
    rhoE_bc = rho_bc*(_eos.e_from_p_rho(press_bc, rho_bc)+0.5*vel_bc*vel_bc);
    Real temp = _eos.temperature_from_p_rho(press_bc, rho_bc);
    epsilon_bc = _ics.a()*temp*temp*temp*temp;
//    }
  }

  // Switch statement for type of equation:
  Real eps = 0.;
  switch (_eqn_type)
  {
  case continuity:
    return rho_bc*vel_bc*_normals[_qp](0)*_test[_i][_qp];
    break;
  case x_momentum:
    return (rhou_bc*vel_bc +press_bc + epsilon_bc/3)*_normals[_qp](0)*_test[_i][_qp];
    break;
  case energy:
    return (rhoE_bc+press_bc)*vel_bc*_normals[_qp](0)*_test[_i][_qp];
    break;
  case radiation:
      eps = epsilon_bc;// -0.5*(_normals[_qp](0)-1)*_ics.eps_hat_pre() + 0.5*(1+_normals[_qp](0))*_ics.eps_hat_post();
//      eps = 0.5*(1+_normals[_qp](0))*_ics.eps_hat_post();
      return (4*vel_bc*epsilon_bc/3*_normals[_qp](0)+0.5*(epsilon_bc-eps))*_test[_i][_qp];
//      return (4*vel_bc*epsilon_bc/3*_normals[_qp](0))*_test[_i][_qp];
//      return 0.5*(_u[_qp]-eps)*_test[_i][_qp];
    break;
  default:
    mooseError("The equation name is not supported in the \"RheaBCs\" type of boundary condition.");
    break;
  }
}

Real
RheaBCs::computeQpJacobian()
{
  return 0;
}

Real
RheaBCs::computeQpOffDiagJacobian(unsigned jvar)
{
  return 0;
}