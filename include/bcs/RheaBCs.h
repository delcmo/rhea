#ifndef RHEABCS_H
#define RHEABCS_H

#include "IntegratedBC.h"
#include "IdealGasEquationOfState.h"
#include "InputFileSpecifiedICsRadHydro.h"

// Forward Declarations
class RheaBCs;
class EquationOfState;

template<>
InputParameters validParams<RheaBCs>();

class RheaBCs : public IntegratedBC
{

public:
  RheaBCs(const InputParameters & parameters);

  virtual ~RheaBCs(){}

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned jvar);

  // Pre-shock parameters
  Real _press_hat_pre;

  // Post-shock parameters
  Real _press_hat_post;

  // Equation name
  enum EFlowEquationType
  {
    continuity = 0,
    x_momentum = 1,
    energy = 2,
    radiation = 3
  };
  MooseEnum _eqn_type;

  // Coupled aux variables
  const VariableValue & _rho;
  const VariableValue & _rho_old;
  const VariableValue & _rhou;
  const VariableValue & _rhou_old;
  const VariableValue & _epsilon;
  const VariableValue & _epsilon_old;
  const VariableValue & _pressure;
  const VariableValue & _pressure_old;

  // Material property:
  const MaterialProperty<Real> & _D;

  // Equation of state
  const IdealGasEquationOfState & _eos;

  // Userobject computing the ICs
  const InputFileSpecifiedICsRadHydro & _ics;

  // Non-dimensional numbers
  const Real _Po;
  const Real _K;

  // Integers jacobian terms
  unsigned int _rho_nb;
  unsigned int _rhou_nb;
  unsigned int _rhoE_nb;
  unsigned int _epsilon_nb;
};

#endif // RHEABCS_H