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

#ifndef INTERNALENERGYAUX_H
#define INTERNALENERGYAUX_H

#include "AuxKernel.h"

class InternalEnergyAux;

template<>
InputParameters validParams<InternalEnergyAux>();

class InternalEnergyAux : public AuxKernel
{
public:

  InternalEnergyAux(const InputParameters & parameters);

protected:
  virtual Real computeValue();

  // Coupled variables:
  const VariableValue & _rho;
  const VariableValue & _rhou;
  const VariableValue & _rhoE;
};

#endif //INTERNALENERGYAUX_H