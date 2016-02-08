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

#ifndef RHEAIC_H
#define RHEAIC_H

// MOOSE Includes
#include "InitialCondition.h"
#include "IdealGasEquationOfState.h"
#include "InputFileSpecifiedICsRadHydro.h"

// Forward Declarations
class RheaIC;

template<>
InputParameters validParams<RheaIC>();

/**
 * RheaIC just returns a constant value.
 */
class RheaIC : public InitialCondition
{
public:

  /**
   * Constructor: Same as the rest of the MOOSE Objects
   */
  RheaIC(const InputParameters & parameters);

  virtual Real value(const Point & p);

private:

  // Position of the membrane:
  Real _membrane;
  Real _length;

  // Equation of state:
  const IdealGasEquationOfState & _eos;

  // Userobject computing the ICs
  const InputFileSpecifiedICsRadHydro & _ics;
};

#endif // RHEAIC_H