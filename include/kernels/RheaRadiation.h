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

class RheaRadiation;

template<>
InputParameters validParams<RheaRadiation>();
class RheaRadiation : public Kernel
{
public:

  RheaRadiation(const std::string & name,
             InputParameters parameters);

protected:
 
  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

  virtual Real computeQpOffDiagJacobian( unsigned int jvar );

private:
    // Material coupled variables:
    VariableValue & _vel;
    VariableValue & _temp;
    
    // Material properties:
    MaterialProperty<Real> & _sigma_a;
    MaterialProperty<Real> & _D;
    
    // Constant: speed of light.
    Real _c;
    Real _a;
};

#endif // RHEARADIATION_H
