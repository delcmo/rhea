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

#include "RheaIC.h"

template<>
InputParameters validParams<RheaIC>()
{
  InputParameters params = validParams<InitialCondition>();

  // Membrane position:
  params.addParam<Real>("membrane", 0.5, "The value of the membrane");
  params.addParam<Real>("length", 0., "To smooth the IC over a given length");    
  // Equation of state
  params.addRequiredParam<UserObjectName>("eos", "parameters for eos.");
  // Userobject computing the ICs
  params.addRequiredParam<UserObjectName>("ics", "parameters for ics.");
  
  return params;
}

RheaIC::RheaIC(const InputParameters & parameters) :
    InitialCondition(parameters),
    // Position of the membrane:
    _membrane(getParam<Real>("membrane")),
    _length(getParam<Real>("length")),
    // Equation of state
    _eos(getUserObject<IdealGasEquationOfState>("eos")),
    // Userobject computing the ICs
    _ics(getUserObject<ComputeICsRadHydro>("ics"))
{
  // Check the value of Cv used in the equation of state
  Real a_hat_0 = std::sqrt(_ics.a()*_ics.T_hat_pre()*_ics.T_hat_pre()*_ics.T_hat_pre()*_ics.T_hat_pre()/(_ics.rho_hat_pre()*_ics.P()));
  Real Cv = a_hat_0*a_hat_0/(_ics.T_hat_pre()*_eos.gamma()*(_eos.gamma()-1.));
  if (std::fabs((Cv-_eos.Cv())/Cv)>1.e-6)
  {
    std::cout<<"--------------------------------------------------------------"<<std::endl;
    std::cout<<"--------------------------------------------------------------"<<std::endl;
    mooseWarning("'"<<this->name()<<"': different Cv coefficients computed. Relative error:"<<std::fabs((Cv-_eos.Cv())/Cv)<<". Cv_{computed}="<<Cv<<" and Cv_{eos}="<<_eos.Cv()<<".");
    std::cout<<"--------------------------------------------------------------"<<std::endl;
    std::cout<<"--------------------------------------------------------------"<<std::endl;    
  }
}

Real
RheaIC::value(const Point & p)
{
  // Get the name of the variable this object acts on
  std::string _name_var = _var.name();
      
  // Define and compute parameters used to smooth the initial condition if required
  Real x1 = _membrane - 0.5 * _length;
  Real x2 = x1 + _length;
  Real a_rho, b_rho, a_vel, b_vel, a_t, b_t, a_eps, b_eps;
  a_rho = ( _ics.rho_hat_pre() - _ics.rho_hat_post()) / _length;
  b_rho = ( x1*_ics.rho_hat_post() - x2*_ics.rho_hat_pre() ) / _length;
  a_vel = ( _ics.vel_hat_pre() - _ics.vel_hat_post()) / _length;
  b_vel = ( x1*_ics.vel_hat_post() - x2*_ics.vel_hat_pre() ) / _length;
  a_t = ( _ics.T_hat_pre() - _ics.T_hat_post()) / _length;
  b_t = ( x1*_ics.T_hat_post() - x2*_ics.T_hat_pre() ) / _length;
  a_eps = ( _ics.eps_hat_pre() - _ics.eps_hat_post()) / _length;
  b_eps = ( x1*_ics.eps_hat_post() - x2*_ics.eps_hat_pre()) / _length;

  // Compute the pressure, velocity and temperature values
  Real rho, temp, vel, epsilon;
  if ( p(0) <= x1 )
  {
    rho = _ics.rho_hat_pre();
    vel = _ics.vel_hat_pre();
    temp = _ics.T_hat_pre();
    epsilon = _ics.eps_hat_pre();
  }
  else if ( p(0) > x2 )
  {
    rho = _ics.rho_hat_post();
    vel = _ics.vel_hat_post();
    temp = _ics.T_hat_post();
    epsilon = _ics.eps_hat_post();
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

  // Return the value of the initial condition.
  if ( _name_var == "rho" ) // Density: rho
    return rho;
  else if ( _name_var == "rhou" ) // Momentum: rhou
    return rhou;
  else if ( _name_var == "rhoE" ) // total energy: rhoE
    return rhoE;
  else if ( _name_var == "epsilon" ) // radiation intensity: epsilon
    return epsilon;
  else
    return 0;
}