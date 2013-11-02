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
  	// User Objects
    _eos(getUserObject<EquationOfState>("eos"))
{}

Real
InitialConditions::value(const Point & p)
{
// Get the name of the variable this object acts on
std::string _name_var = _var.name();
    
// Compute the pressure, velocity and temperature values
Real _rho = 0.;
Real _temp = 0.;
Real _vel = 0.;
Real _epsilon = 0.;
if ( p(0) <= _membrane ) {
    _rho = _rho_left;
    _vel = _v_left;
    _temp = _t_left;
    _epsilon = _eps_left;
    }
else if ( p(0) > _membrane ) {
    _rho = _rho_right;
    _vel = _v_right;
    _temp = _t_right;
    _epsilon = _eps_right;
    }
else {
    mooseError("The value of the membrane is probably wrong.");
}
      
// Compute the conservative variables
Real _pressure = _rho*(_eos.Cv()*(_eos.gamma()-1)*_temp) - _eos.Pinf();
Real _int_energy = (_pressure+_eos.gamma()*_eos.Pinf())/(_rho*(_eos.gamma()-1)) + _eos.qcoeff();
Real _tot_energy = _rho*(_int_energy + 0.5*_vel*_vel);
    
// Return the value of the initial condition. Identify the name of the variable
// Density: rhoA
if ( _name_var == "rho" )
	return ( _rho);
// Momentum: rhouA and rhovA
else if ( _name_var == "rhou" )
{
	return ( _rho*_vel);
}
// total energy: rhoE
else if ( _name_var == "rhoE" )
	return ( _tot_energy );
// radiation intensity: epsilon
else if ( _name_var == "epsilon" )
	return ( _epsilon );
else
	return 0;
}
