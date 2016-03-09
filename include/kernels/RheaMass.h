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

#ifndef RHEAMASS_H
#define RHEAMASS_H

#include "Kernel.h"

class RheaMass;

template<>
InputParameters validParams<RheaMass>();
class RheaMass : public Kernel
{
public:

  RheaMass(const InputParameters & parameters);

protected:
 
  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

  virtual Real computeQpOffDiagJacobian( unsigned int jvar );

private:
  // Coupled variables
  const VariableValue & _rhou;

  // Integer for jacobian terms
  unsigned int _rhou_nb;
};

#endif // RheaMass_H
