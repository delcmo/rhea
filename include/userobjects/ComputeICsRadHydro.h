#ifndef COMPUTEICSRADHYDRO_H
#define COMPUTEICSRADHYDRO_H

#include "GeneralUserObject.h"
#include "IdealGasEquationOfState.h"

// Forward Declarations
class ComputeICsRadHydro;

template<>
InputParameters validParams<ComputeICsRadHydro>();

class ComputeICsRadHydro : public GeneralUserObject
{
public:
  // Constructor
  ComputeICsRadHydro(const std::string & name, InputParameters parameters);

  // Destructor
  virtual ~ComputeICsRadHydro();

  /**
   * Called when this object needs to compute something.
   */
  virtual void execute() {}

  /**
   * Called before execute() is ever called so that data can be cleared.
   */
  virtual void initialize(){}


  virtual void destroy();

  virtual void finalize() {};

  // Pre- and post-shock functions
  Real rho_hat_post() const { return _rho_hat_post;};
  Real rho_hat_pre() const { return _rho_hat_pre;};
  Real T_hat_post() const { return _T_hat_post;};
  Real T_hat_pre() const { return _T_hat_pre;};
  Real vel_hat_post() const { return _vel_hat_post;};
  Real vel_hat_pre() const { return _vel_hat_pre;};
  Real eps_hat_post() const { return _eps_hat_post;};
  Real eps_hat_pre() const { return _eps_hat_pre;};

  // Non-dimensionalized functions
  Real P() const { return _P;};
  Real K() const { return _K;};
  Real SIGMA_A() const { return _SIGMA_A;};
  Real a() const { return _a;};
  Real c() const { return _sp;};

protected:
  // Physical properties
  Real _sp;
  Real _a;

  // Non-dimensionalized numbers
  Real _P;
  Real _Mach_inlet;
  Real _K;
  Real _SIGMA_A;

  // Pre-shock parameters
  Real _rho_hat_pre;
  Real _T_hat_pre;
  Real _vel_hat_pre;
  Real _eps_hat_pre;

  // Post-shock parameters
  Real _rho_hat_post;
  Real _T_hat_post;
  Real _vel_hat_post;
  Real _eps_hat_post;

  // Equation of state
  const IdealGasEquationOfState & _eos;
};


#endif // COMPUTEICSRADHYDRO_H