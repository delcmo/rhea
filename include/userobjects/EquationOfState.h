#ifndef EQUATIONOFSTATE_H
#define EQUATIONOFSTATE_H

#include "GeneralUserObject.h"

// Forward Declarations
class EquationOfState;

template<>
InputParameters validParams<EquationOfState>();

class EquationOfState : public GeneralUserObject
{
public:
  // Constructor
  EquationOfState(const std::string & name, InputParameters parameters);

  // Destructor
  virtual ~EquationOfState();

  /**
   * Called when this object needs to compute something.
   */
  virtual void execute() {};

  /**
   * Called before execute() is ever called so that data can be cleared.
   */
  virtual void initialize(){};

  /**
   * Finalize.  This is called _after_ execute() and _after_ threadJoin()!  This is probably where you want to do MPI communication!
   */
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

  // The derivative of pressure wrt density (rho)
  virtual Real dp_drho(Real rho, Real rhou, Real rhoE) const;

  // The derivative of pressure wrt x-momentum (rhou)
  virtual Real dp_drhou(Real rho, Real vel, Real rhoE) const;

  // The derivative of pressure wrt energy (rhoE)
  virtual Real dp_drhoE(Real rho, Real vel, Real rhoE) const;

  // The derivative of temperature wrt density (rho)
  virtual Real dT_drho(Real rho, Real rhou, Real rhoE) const;

  // The derivative of temperature wrt x-momentum (rhou)
  virtual Real dT_drhou(Real rho, Real vel, Real rhoE) const;

  // The derivative of temperature wrt energy (rhoE)
  virtual Real dT_drhoE(Real rho, Real vel, Real rhoE) const;

  // Sound speed squared
  virtual Real c2(Real rho, Real vel, Real rhoE, Real epsilon) const;

  // Sound speed squared
  virtual Real c2_from_p_rho(Real rho, Real pressure, Real epsilon) const;

protected:
  // Prints an error message for non-implemented functions
  void error_not_implemented(std::string method_name) const;
};


#endif // EQUATIONOFSTATE_H

