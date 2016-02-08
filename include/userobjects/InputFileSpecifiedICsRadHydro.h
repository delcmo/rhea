#ifndef INPUTFILESPECIFIEDICSRADHYDRO_H
#define INPUTFILESPECIFIEDICSRADHYDRO_H

#include "GeneralUserObject.h"
#include "IdealGasEquationOfState.h"

// Forward Declarations
class InputFileSpecifiedICsRadHydro;

template<>
InputParameters validParams<InputFileSpecifiedICsRadHydro>();

class InputFileSpecifiedICsRadHydro : public GeneralUserObject
{
public:
  // Constructor
  InputFileSpecifiedICsRadHydro(const InputParameters & parameters);

  // Destructor
  virtual ~InputFileSpecifiedICsRadHydro();

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

  // Boolean for dimensional form
  bool _is_dmsl_form;

  // Pre- and post-shock functions
  Real rho_hat_post() const { return _rho_hat_post;};
  Real rho_hat_pre() const { return _rho_hat_pre;};
  Real T_hat_post() const { return _T_hat_post;};
  Real T_hat_pre() const { return _T_hat_pre;};
  Real vel_hat_post() const { return _vel_hat_post;};
  Real vel_hat_pre() const { return _vel_hat_pre;};
  Real eps_hat_post() const { return _eps_hat_post;};
  Real eps_hat_pre() const { return _eps_hat_pre;};
  Real mach_hat_post() const { return _mach_hat_post;};
  Real mach_hat_pre() const { return _mach_hat_pre;};

  // Non-dimensionalized functions
  Real P() const { return _P;};
  Real K() const { return _K;};
  Real SIGMA_A() const { return _SIGMA_A;};
  Real C() const { return _C;};
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
  Real _C;

  // Pre-shock parameters
  Real _rho_hat_pre;
  Real _T_hat_pre;
  Real _vel_hat_pre;
  Real _eps_hat_pre;
  Real _mach_hat_pre;

  // Post-shock parameters
  Real _rho_hat_post;
  Real _T_hat_post;
  Real _vel_hat_post;
  Real _eps_hat_post;
  Real _mach_hat_post;

  // Equation of state
  const IdealGasEquationOfState & _eos;
};


#endif // INPUTFILESPECIFIEDICSRADHYDRO_H
