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

#ifndef MACHNUMBERAUX_H
#define MACHNUMBERAUX_H

#include "AuxKernel.h"
#include "EquationOfState.h"


//Forward Declarations
class MachNumberAux;

template<>
InputParameters validParams<MachNumberAux>();

/**
 * Coupled auxiliary value
 */
class MachNumberAux : public AuxKernel
{
public:

  MachNumberAux(const InputParameters & parameters);

protected:
  virtual Real computeValue();

  // Coupled variables
  const VariableValue & _rho;
  const VariableValue & _rhou;
  const VariableValue & _rhoE;
  const VariableValue & _epsilon;
  // Equation of state
  const EquationOfState & _eos;
};

#endif //MachNumberAux_H