#ifndef ENTROPYVISCOSITYCOEFFICIENT_H
#define ENTROPYVISCOSITYCOEFFICIENT_H

#include "Material.h"
#include "MaterialProperty.h"
#include "EquationOfState.h"

//Forward Declarations
class EntropyViscosityCoefficient;

template<>
InputParameters validParams<EntropyViscosityCoefficient>();

class EntropyViscosityCoefficient : public Material
{
public:
  EntropyViscosityCoefficient(const std::string & name, InputParameters parameters);

protected:
  virtual void computeQpProperties();

private:

  // Coupled variables
  VariableValue & _rho;
  VariableValue & _rho_old;
  VariableValue & _rho_older;
  VariableGradient & _grad_rho;
  VariableValue & _rhou;
  VariableValue & _eps;
  VariableValue & _eps_old;
  VariableValue & _eps_older;
  VariableGradient & _grad_eps;

  // Coupled aux variables
  VariableValue & _press;
  VariableValue & _press_old;
  VariableValue & _press_older;
  VariableGradient & _grad_press;

  // Jump values:
  VariableValue & _jump_press;
  VariableValue & _jump_dens;

  // Material properties: viscosity coefficients.
  MaterialProperty<Real> & _kappa;
  MaterialProperty<Real> & _kappa_max;

  // Multiplicative coefficient for viscosity:
  double _Cjump;
  bool _is_first_order_viscosity;

  // UserObject: equation of state
  const EquationOfState & _eos;
};

#endif // ENTROPYVISCOSITYCOEFFICIENT_H