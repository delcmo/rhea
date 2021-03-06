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

#ifndef RHEAENERGY_H
#define RHEAENERGY_H

#include "Kernel.h"
#include "EquationOfState.h"
#include "InputFileSpecifiedICsRadHydro.h"

class RheaEnergy;

template<>
InputParameters validParams<RheaEnergy>();
class RheaEnergy : public Kernel
{
public:

  RheaEnergy(const InputParameters & parameters);

protected:
 
  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

  virtual Real computeQpOffDiagJacobian( unsigned int jvar );

private:

  // Booleans for dimensional form
  bool _is_dmsl_form;

  // coupled variables:
  const VariableValue & _rho;
  const VariableValue & _rhou;
  const VariableValue & _epsilon;
  const VariableGradient & _grad_eps;

  // Equation of state
  const EquationOfState & _eos;

  // Material property:
  const MaterialProperty<Real> & _sigma_a;

  // Userobject computing the ICs
  const InputFileSpecifiedICsRadHydro & _ics;

  // Dimensional numbers
  const Real _c;
  const Real _a;

  // Non-dimensional numbers
  const Real _Po;
  const Real _SIGMA;

  // Integers jacobian terms
  unsigned int _rho_nb;
  unsigned int _rhou_nb;
  unsigned int _epsilon_nb;
};

#endif // RHEAENERGY_H
