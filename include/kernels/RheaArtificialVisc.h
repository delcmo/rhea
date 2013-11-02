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
        CONTINUITY = 0,
        MOMENTUM = 1,
        ENERGY = 2,
        RADIATION = 3
    };

    // Diffusion name
    std::string _equ_name;
    
    // Diffusion type
    MooseEnum _equ_type;
    
    // Boolean for dissipation in radiation:
    bool _isRadiation;
    
    // Material variables:
    VariableValue & _rho;
    VariableGradient & _grad_rho;
    VariableValue & _vel;
    VariableGradient & _grad_vel;
    VariableGradient & _grad_rhoe;
    
    // Radiation variables:
    VariableGradient & _grad_epsilon;
    
    // Material property: viscosity coefficient.
    MaterialProperty<Real> & _mu;
    MaterialProperty<Real> & _kappa;
};

#endif //RHEAARTIFICIALVISC_H
