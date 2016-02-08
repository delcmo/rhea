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

#ifndef RHEANONDIMENSIONALENERGY_H
#define RHEANONDIMENSIONALENERGY_H

#include "Kernel.h"
#include "EquationOfState.h"
#include "InputFileSpecifiedICsRadHydro.h"

class RheaNonDimensionalEnergy;

template<>
InputParameters validParams<RheaNonDimensionalEnergy>();
class RheaNonDimensionalEnergy : public Kernel
{
public:

  RheaNonDimensionalEnergy(const InputParameters & parameters);

protected:
 
  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

  virtual Real computeQpOffDiagJacobian( unsigned int jvar );

private:

  // coupled variables:
  VariableValue & _rho;
  VariableValue & _rhou;
  VariableValue & _epsilon;
  VariableGradient & _grad_eps;

  // Equation of state
  const EquationOfState & _eos;

  // Material property:
  const MaterialProperty<Real> & _sigma_a;

  // Userobject computing the ICs
  const InputFileSpecifiedICsRadHydro & _ics;

  // Non-dimensional numbers
  Real _Po;
  Real _SIGMA;

  // Integers jacobian terms
  unsigned int _rho_nb;
  unsigned int _rhou_nb;
  unsigned int _epsilon_nb;
};

#endif // RHEANONDIMENSIONALENERGY_H
