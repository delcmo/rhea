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

#ifndef RHEANONDIMENSIONALRADIATION_H
#define RHEANONDIMENSIONALRADIATION_H

#include "Kernel.h"
#include "EquationOfState.h"
#include "InputFileSpecifiedICsRadHydro.h"

class RheaNonDimensionalRadiation;

template<>
InputParameters validParams<RheaNonDimensionalRadiation>();
class RheaNonDimensionalRadiation : public Kernel
{
public:

  RheaNonDimensionalRadiation(const InputParameters & parameters);

protected:
 
  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

  virtual Real computeQpOffDiagJacobian( unsigned int jvar );

private:

  // Coupled variables:
  const VariableValue & _rho;
  const VariableValue & _rhou;
  const VariableValue & _rhoE;

  // Equation of state
  const EquationOfState & _eos;

  // Material properties:
  const MaterialProperty<Real> & _sigma_a;
  const MaterialProperty<Real> & _D;

  // Userobject computing the ICs
  const InputFileSpecifiedICsRadHydro & _ics;

  // Non-dimensional numbers
  const Real _SIGMA;
  const Real _K;

  // Integers for jacobian terms
  unsigned int _rho_nb;
  unsigned int _rhou_nb;
  unsigned int _rhoE_nb;
};

#endif // RHEANONDIMENSIONALRADIATION_H