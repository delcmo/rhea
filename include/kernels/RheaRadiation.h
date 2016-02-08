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

#ifndef RHEARADIATION_H
#define RHEARADIATION_H

#include "Kernel.h"
#include "EquationOfState.h"
#include "InputFileSpecifiedICsRadHydro.h"

class RheaRadiation;

template<>
InputParameters validParams<RheaRadiation>();
class RheaRadiation : public Kernel
{
public:

  RheaRadiation(const InputParameters & parameters);

protected:
 
  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

  virtual Real computeQpOffDiagJacobian( unsigned int jvar );

private:

  // Booleans for dimensional form
  bool _is_dmsl_form;

  // Coupled variables:
  VariableValue & _rho;
  VariableValue & _rhou;
  VariableValue & _rhoE;

  // Equation of state
  const EquationOfState & _eos;

  // Material properties:
  const MaterialProperty<Real> & _sigma_a;
  const MaterialProperty<Real> & _D;

  // Userobject computing the ICs
  const InputFileSpecifiedICsRadHydro & _ics;

  // Dimensional numbers
  Real _c;
  Real _a;

  // Non-dimensional numbers
  Real _SIGMA;
  Real _K;

  // Integers for jacobian terms
  unsigned int _rho_nb;
  unsigned int _rhou_nb;
  unsigned int _rhoE_nb;
};

#endif // RHEARADIATION_H
