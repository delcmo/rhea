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

#ifndef RADTEMPAUX_H
#define RADTEMPAUX_H

#include "AuxKernel.h"
#include "EquationOfState.h"

//Forward Declarations
class RadTempAux;

template<>
InputParameters validParams<RadTempAux>();

class RadTempAux : public AuxKernel
{
public:

  RadTempAux(const InputParameters & parameters);

protected:
  virtual Real computeValue();

  // Coupled variables
  const VariableValue & _eps;
  // Constant
  const Real _a;
};

#endif //RADTEMPAUX_H