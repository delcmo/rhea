#include "IdealGasEquationOfState.h"
#include "MooseError.h"

template<>
InputParameters validParams<IdealGasEquationOfState>()
{
  InputParameters params = validParams<EquationOfState>();

  params.addParam<Real>("gamma", 1.4, "  gamma");
  params.addParam<Real>("Cv", 1., "Heat coefficient");
  
  return params;
}

IdealGasEquationOfState::IdealGasEquationOfState(const InputParameters & parameters) :
  EquationOfState(parameters),
    _gamma(getParam<Real>("gamma")),
    _Cv(getParam<Real>("Cv"))
{
}

IdealGasEquationOfState::~IdealGasEquationOfState()
{
}

void
IdealGasEquationOfState::destroy()
{
}

Real IdealGasEquationOfState::pressure(Real rho, Real vel, Real rhoE) const
{
  Real rhoe = rhoE - 0.5 * rho*vel*vel;
  return (_gamma-1) * rhoe;
}

Real IdealGasEquationOfState::temperature(Real rho, Real vel, Real rhoE) const
{
  Real pressure = this->pressure(rho, vel, rhoE);
  return this->temperature_from_p_rho(pressure, rho);
}

Real IdealGasEquationOfState::rho_from_p_T(Real pressure, Real temperature) const
{
  return pressure / ((_gamma-1)*_Cv*temperature);
}

Real IdealGasEquationOfState::e_from_p_rho(Real pressure, Real rho) const
{
  return pressure / ((_gamma-1)*rho);
}

Real IdealGasEquationOfState::temperature_from_p_rho(Real pressure, Real rho) const
{
  return pressure / ((_gamma-1)*_Cv*rho);
}

Real IdealGasEquationOfState::p_from_T_rho(Real temperature, Real rho) const
{
  return temperature*(_gamma-1)*_Cv*rho;
}

Real IdealGasEquationOfState::dp_drho(Real rho, Real rhou, Real rhoE) const
{
  return 0.;
}

Real IdealGasEquationOfState::dp_drhou(Real rho, Real rhou, Real rhoE) const
{
  return 0.;
}

Real IdealGasEquationOfState::dp_drhoE(Real rho, Real rhou, Real rhoE) const
{
  return 0.;
}

Real IdealGasEquationOfState::dT_drho(Real rho, Real rhou, Real rhoE) const
{
  return 0.;
}

Real IdealGasEquationOfState::dT_drhou(Real rho, Real rhou, Real rhoE) const
{
  return 0.;
}

Real IdealGasEquationOfState::dT_drhoE(Real rho, Real rhou, Real rhoE) const
{
  return 0.;
}

Real IdealGasEquationOfState::c2(Real rho, Real vel, Real rhoE, Real epsilon) const
{
  Real pressure = this->pressure(rho, vel, rhoE);
  return this->c2_from_p_rho(rho, pressure, epsilon);
}

Real IdealGasEquationOfState::c2_from_p_rho(Real rho, Real pressure, Real epsilon) const
{
  return (_gamma * pressure +4.*epsilon/9.) / rho;
}