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

#include "InitialConditions.h"

template<>
InputParameters validParams<InitialConditions>()
{
  InputParameters params = validParams<InitialCondition>();

  // Initial conditions:
  params.addRequiredParam<Real>("rho_init_left", "Initial density on the left");
  params.addRequiredParam<Real>("rho_init_right", "Initial density on the right");
  params.addRequiredParam<Real>("vel_init_left", "Initial velocity on the left");
  params.addRequiredParam<Real>("vel_init_right", "Inital velocity on the right");
  params.addRequiredParam<Real>("temp_init_left", "Initial value of the temperature");
  params.addRequiredParam<Real>("temp_init_right", "Initial value of the temperature");
  params.addRequiredParam<Real>("eps_init_left", "Initial value of the radiation intensity.");
  params.addRequiredParam<Real>("eps_init_right", "Initial value of the radiation intensity");
  // Membrane position:
  params.addParam<Real>("membrane", 0.5, "The value of the membrane");
  params.addParam<Real>("length", 0., "To smooth the IC over a given length");    
  // Equation of state
  params.addRequiredParam<UserObjectName>("eos", "parameters for eos.");
  
  return params;
}

InitialConditions::InitialConditions(const std::string & name,
                     InputParameters parameters) :
    InitialCondition(name, parameters),
    // IC parameters
    _rho_left(getParam<Real>("rho_init_left")),
    _rho_right(getParam<Real>("rho_init_right")),
    _v_left(getParam<Real>("vel_init_left")),
    _v_right(getParam<Real>("vel_init_right")),
    _t_left(getParam<Real>("temp_init_left")),
    _t_right(getParam<Real>("temp_init_right")),
    _eps_left(getParam<Real>("eps_init_left")),
    _eps_right(getParam<Real>("eps_init_right")),
    // Position of the membrane:
    _membrane(getParam<Real>("membrane")),
    _length(getParam<Real>("length")),
    // User Objects
    _eos(getUserObject<EquationOfState>("eos"))
{}

Real
InitialConditions::value(const Point & p)
{
  // Get the name of the variable this object acts on
  std::string _name_var = _var.name();
      
  // Define and compute parameters used to smooth the initial condition if wished
  Real x1 = _membrane - 0.5 * _length;
  Real x2 = x1 + _length;
  Real a_rho, b_rho, a_vel, b_vel, a_t, b_t, a_eps, b_eps;
  a_rho = ( _rho_left - _rho_right) / _length;
  b_rho = ( x1*_rho_right - x2*_rho_left ) / _length;
  a_vel = ( _v_left - _v_right) / _length;
  b_vel = ( x1*_v_right - x2*_v_left ) / _length;
  a_t = ( _t_left - _t_right) / _length;
  b_t = ( x1*_t_right - x2*_t_left ) / _length;
  a_eps = ( _eps_left - _eps_right) / _length;
  b_eps = ( x1*_eps_right - x2*_eps_left ) / _length;

  // Compute the pressure, velocity and temperature values
  Real rho = 0.;
  Real temp = 0.;
  Real vel = 0.;
  Real epsilon = 0.;
  if ( p(0) <= x1 )
  {
    rho = _rho_left;
    vel = _v_left;
    temp = _t_left;
    epsilon = _eps_left;
  }
  else if ( p(0) > x2 )
  {
    rho = _rho_right;
    vel = _v_right;
    temp = _t_right;
    epsilon = _eps_right;
  }
  else
  {
    rho = ( a_rho * p(0) + b_rho );
    vel = ( a_vel * p(0) + b_vel );
    temp = ( a_t * p(0) + b_t );
    epsilon = ( a_eps * p(0) + b_eps );
  }

  // Compute the conservative variables
  Real rhou = rho*vel;
  Real press = _eos.p_from_T_rho(temp, rho);
  Real e = _eos.e_from_p_rho(press, rho);
  Real rhoE = rho*(e + 0.5*vel*vel);

  // Return the value of the initial condition. Identify the name of the variable
  // Density: rho
  if ( _name_var == "rho" )
    return rho;
  // Momentum: rhou
  else if ( _name_var == "rhou" )
    return rhou;
  // total energy: rhoE
  else if ( _name_var == "rhoE" )
    return rhoE;
  // radiation intensity: epsilon
  else if ( _name_var == "epsilon" )
    return epsilon;
  else
    return 0;
}