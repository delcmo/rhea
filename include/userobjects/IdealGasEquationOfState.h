#ifndef IDEALGASEQUATIONOFSTATE_H
#define IDEALGASEQUATIONOFSTATE_H

#include "EquationOfState.h"

class IdealGasEquationOfState;

template<>
InputParameters validParams<IdealGasEquationOfState>();

class IdealGasEquationOfState : public EquationOfState
{
public:
  // Constructor
  IdealGasEquationOfState(const std::string & name, InputParameters parameters);

  // Destructor  
  virtual ~IdealGasEquationOfState(); 

  virtual void execute() {};

  virtual void initialize() {};


  virtual void destroy();

  virtual void finalize() {};

  // The interface for derived EquationOfState objects to implement...
  virtual Real pressure(Real rho, Real vel, Real rhoE) const;

  // The interface for derived EquationOfState objects to implement...
  virtual Real temperature(Real rho, Real vel, Real rhoE) const;

  // density from pressure and temperature
  virtual Real rho_from_p_T(Real pressure, Real temperature) const;

  // internal energy from pressure and density
  virtual Real e_from_p_rho(Real pressure, Real rho) const;

  // temperature from pressure and density
  virtual Real temperature_from_p_rho(Real pressure, Real rho) const;

  // pressure from temperature and density
  virtual Real p_from_T_rho(Real temperature, Real rho) const;

  // The derivative of pressure wrt density (h)
  virtual Real dp_drho(Real rho, Real rhou, Real rhoE) const;

  // The derivative of pressure wrt x-momentum (rhou)
  virtual Real dp_drhou(Real rho, Real vel, Real rhoE) const;

  // Sound speed squared
  virtual Real c2(Real rho, Real vel, Real rhoE, Real epsilon) const;

  // Sound speed squared
  virtual Real c2_from_p_rho(Real rho, Real pressure, Real epsilon) const;

  Real gamma() const { return _gamma; };

  Real Cv() const { return _Cv; };

protected:
  Real _gamma;
  Real _Cv;
};

#endif // IDEALGASEQUATIONOFSTATE_H