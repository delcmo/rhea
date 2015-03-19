#include "EquationOfState.h"
#include "MooseError.h"

template<>
InputParameters validParams<EquationOfState>()
{
  InputParameters params = validParams<UserObject>();
  return params;
}

EquationOfState::EquationOfState(const std::string & name, InputParameters parameters) :
    GeneralUserObject(name, parameters)
{}

EquationOfState::~EquationOfState()
{
}

// The pressure as a function of the conservative variables
Real EquationOfState::pressure(Real, Real, Real) const
{
  this->error_not_implemented("pressure");
  return 0.;
}

// The temperature as a function of the conservative variables
Real EquationOfState::temperature(Real, Real, Real) const
{
  this->error_not_implemented("temperature");
  return 0.;
}

// The density as a function of pressure and temperature
Real EquationOfState::rho_from_p_T(Real, Real) const
{
  this->error_not_implemented("rho from pressure and temperature");
  return 0.;
}

// The temperature as a function of pressure and density
Real EquationOfState::e_from_p_rho(Real, Real) const
{
  this->error_not_implemented("temperature");
  return 0.;
}

// The temperature as a function of pressure and density
Real EquationOfState::temperature_from_p_rho(Real, Real) const
{
  this->error_not_implemented("temperature");
  return 0.;
}

// pressure from temperature and density
Real EquationOfState::p_from_T_rho(Real temperature, Real rho) const
{
  this->error_not_implemented("p_from_T_rho");
  return 0.;
}

// The derivative of pressure wrt rho
Real EquationOfState::dp_drho(Real, Real, Real) const
{
  this->error_not_implemented("dp_drho");
  return 0.;
}

// The derivative of pressure wrt rhou
Real EquationOfState::dp_drhou(Real, Real, Real) const
{
  this->error_not_implemented("dp_drhou");
  return 0.;
}

// Sound speed square
Real
EquationOfState::c2(Real, Real, Real, Real) const
{
  this->error_not_implemented("c2");
  return 0.;
}

// Sound speed square
Real
EquationOfState::c2_from_p_rho(Real, Real, Real) const
{
  this->error_not_implemented("c2");
  return 0.;
}

void EquationOfState::error_not_implemented(std::string method_name) const
{
  mooseError("Your EquationOfState object does not implement: " + method_name);
}