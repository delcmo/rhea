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

#ifndef RHEANONDIMENSIONALMOMENTUM_H
#define RHEANONDIMENSIONALMOMENTUM_H

#include "Kernel.h"
#include "EquationOfState.h"
#include "InputFileSpecifiedICsRadHydro.h"

class RheaNonDimensionalMomentum;

template<>
InputParameters validParams<RheaNonDimensionalMomentum>();
class RheaNonDimensionalMomentum : public Kernel
{
public:

  RheaNonDimensionalMomentum(const InputParameters & parameters);

protected:
 
  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

  virtual Real computeQpOffDiagJacobian( unsigned int jvar );

private:
  // coupled variables:
  const VariableValue & _rho;
  const VariableValue & _rhoE;
  const VariableValue & _epsilon;

  // userobjects: EoS and ICs
  const EquationOfState & _eos;
  const InputFileSpecifiedICsRadHydro & _ics;

  // Non-dimensional number Po
  const Real _Po;

  // Integers for jacobian terms
  unsigned int _rho_nb;
  unsigned int _rhoE_nb;
  unsigned int _epsilon_nb;
};

#endif // RHEANONDIMENSIONALMOMENTUM_H