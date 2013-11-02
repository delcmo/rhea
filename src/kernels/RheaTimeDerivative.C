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

#include "RheaTimeDerivative.h"

template<>
InputParameters validParams<RheaTimeDerivative>()
{
  InputParameters params = validParams<TimeDerivative>();
  return params;
}

RheaTimeDerivative::RheaTimeDerivative(const std::string & name,
                                             InputParameters parameters) :
    TimeDerivative(name,parameters)
{}

Real
RheaTimeDerivative::computeQpResidual()
{
  return TimeDerivative::computeQpResidual();
}

Real
RheaTimeDerivative::computeQpJacobian()
{
  return TimeDerivative::computeQpJacobian();
}
