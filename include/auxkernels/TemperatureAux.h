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

#ifndef TEMPERATUREAUX_H
#define TEMPERATUREAUX_H

#include "AuxKernel.h"
#include "EquationOfState.h"


//Forward Declarations
class TemperatureAux;

template<>
InputParameters validParams<TemperatureAux>();

class TemperatureAux : public AuxKernel
{
public:

  TemperatureAux(const InputParameters & parameters);

protected:
  virtual Real computeValue();

  // Coupled variables
  const VariableValue & _rho;
  const VariableValue & _rhou;
  const VariableValue & _rhoE;
  // Equation of state
  const EquationOfState & _eos;
};

#endif //TEMPERATUREAUX_H