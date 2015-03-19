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

#ifndef RHEAARTIFICIALVISC_H
#define RHEAARTIFICIALVISC_H

#include "Kernel.h"

// Forward Declarations
class RheaArtificialVisc;

template<>
InputParameters validParams<RheaArtificialVisc>();

class RheaArtificialVisc : public Kernel
{
public:

  RheaArtificialVisc(const std::string & name,
             InputParameters parameters);

protected:

  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

  virtual Real computeQpOffDiagJacobian(unsigned int _jvar);
    
private:
  // Equations types
  enum EquationType
  {
      continuity = 0,
      x_momentum = 1,
      energy = 2,
      radiation = 3
  };
  MooseEnum _equ_type;

  // Boolean for dissipation in radiation:
  bool _isRadiation;

  // Material property:
  MaterialProperty<Real> & _kappa;
  MaterialProperty<Real> & _D;
};

#endif //RHEAARTIFICIALVISC_H
