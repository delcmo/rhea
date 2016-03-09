#ifndef TIMESTEPCFL_H
#define TIMESTEPCFL_H

#include "ElementPostprocessor.h"
#include "EquationOfState.h"

class TimeStepCFL;

template<>
InputParameters validParams<TimeStepCFL>();

/**
 * The inviscid time step stability limit:
 *
 * h_e \over {|\vec u| + c}
 */
class TimeStepCFL : public ElementPostprocessor
{
public:
  TimeStepCFL(const InputParameters & parameters);
  virtual ~TimeStepCFL();

  virtual void initialize();
  virtual void execute();
  virtual Real getValue();
  virtual void threadJoin(const UserObject & uo);

protected:

  // Coupled variables
  const VariableValue & _rho;
  const VariableValue & _rhou;
  const VariableValue & _rhoE;
  const VariableValue & _epsilon;

  // Equation of state
  const EquationOfState & _eos;

  // Parameter
  const Real _cfl;
  Real _value;
};


#endif // TIMESTEPCFL_H
