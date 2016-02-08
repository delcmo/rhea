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
#include "EquationOfState.h"
#include "InputFileSpecifiedICsRadHydro.h"

class RheaMomentum;

template<>
InputParameters validParams<RheaMomentum>();
class RheaMomentum : public Kernel
{
public:

  RheaMomentum(const InputParameters & parameters);

protected:
 
  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

  virtual Real computeQpOffDiagJacobian( unsigned int jvar );

private:
  // Booleans for dimensional form
  bool _is_dmsl_form;

  // coupled variables:
  VariableValue & _rho;
  VariableValue & _rhoE;
  VariableValue & _epsilon;

  // userobjects: EoS and ICs
  const EquationOfState & _eos;
  const InputFileSpecifiedICsRadHydro & _ics;

  // Non-dimensional number Po
  Real _Po;

  // Integers for jacobian terms
  unsigned int _rho_nb;
  unsigned int _rhoE_nb;
  unsigned int _epsilon_nb;
};

#endif // RHEAMOMENTUM_H