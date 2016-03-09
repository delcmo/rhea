#ifndef ENTROPYVISCOSITYCOEFFICIENT_H
#define ENTROPYVISCOSITYCOEFFICIENT_H

#include "Material.h"
#include "MaterialProperty.h"
#include "EquationOfState.h"
#include "InputFileSpecifiedICsRadHydro.h"

//Forward Declarations
class EntropyViscosityCoefficient;

template<>
InputParameters validParams<EntropyViscosityCoefficient>();

class EntropyViscosityCoefficient : public Material
{
public:
  EntropyViscosityCoefficient(const InputParameters & parameters);

protected:
  virtual void computeQpProperties();

private:
  // Boolean
  bool _is_dmsl_form;

  // Coupled variables
  const VariableValue & _rho;
  const VariableValue & _rho_old;
  const VariableValue & _rho_older;
  const VariableGradient & _grad_rho;
  const VariableValue & _rhou;
  const VariableValue & _eps;
  const VariableValue & _eps_old;
  const VariableValue & _eps_older;
  const VariableGradient & _grad_eps;

  // Coupled aux variables
  const VariableValue & _press;
  const VariableValue & _press_old;
  const VariableValue & _press_older;
  const VariableGradient & _grad_press;

  // Jump values:
  const VariableValue & _jump_press;
  const VariableValue & _jump_dens;

  // Material properties: viscosity coefficients.
  MaterialProperty<Real> & _kappa;
  MaterialProperty<Real> & _kappa_max;

  // Multiplicative coefficient for viscosity:
  const Real _Cjump;
  bool _is_first_order_viscosity;
  bool _use_jumps;

  // UserObject: equation of state
  const EquationOfState & _eos;

  // Userobject computing the ICs
  const InputFileSpecifiedICsRadHydro & _ics;

  // Non-dimensional number
  const Real _Po;
};

#endif // ENTROPYVISCOSITYCOEFFICIENT_H