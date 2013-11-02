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

#ifndef RHEAENERGY_H
#define RHEAENERGY_H

#include "Kernel.h"

class RheaEnergy;

template<>
InputParameters validParams<RheaEnergy>();
class RheaEnergy : public Kernel
{
public:

  RheaEnergy(const std::string & name,
             InputParameters parameters);

protected:
 
  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

  virtual Real computeQpOffDiagJacobian( unsigned int jvar );

private:
    // Material coupled variables:
    VariableValue & _vel;
    VariableValue & _temp;
    VariableValue & _pressure;
    
    // Radiation variable:
    VariableValue & _epsilon;
    VariableGradient & _grad_eps;
    
    // Material property:
    MaterialProperty<Real> & _sigma_a;
    
    // Constant: speed of light.
    Real _c;
    Real _a;
};

#endif // RHEAENERGY_H
