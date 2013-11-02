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

#ifndef RHEAMOMENTUM_H
#define RHEAMOMENTUM_H

#include "Kernel.h"

class RheaMomentum;

template<>
InputParameters validParams<RheaMomentum>();
class RheaMomentum : public Kernel
{
public:

  RheaMomentum(const std::string & name,
             InputParameters parameters);

protected:
 
  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

  virtual Real computeQpOffDiagJacobian( unsigned int jvar );

private:
    // Material coupled variables:
    VariableValue & _vel;
    VariableValue & _pressure;
    
    // Radiation variable:
    VariableValue & _epsilon;
};

#endif // RHEAMOMENTUM_H
