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
/**
This function computes the fluid internal energy 'rhoe' from the conservative variables. It is dimension agnostic.
**/
#include "InternalEnergyAux.h"

template<>
InputParameters validParams<InternalEnergyAux>()
{
  InputParameters params = validParams<AuxKernel>();
    // Coupled variables
    params.addRequiredCoupledVar("rho", "density");
    params.addRequiredCoupledVar("rhou", "momentum");
    params.addRequiredCoupledVar("rhoE", "material total energy");
  return params;
}

InternalEnergyAux::InternalEnergyAux(const std::string & name, InputParameters parameters) :
    AuxKernel(name, parameters),
    // Coupled variables:
    _rho(coupledValue("rho")),
    _rhou(coupledValue("rhou")),
    _rhoE(coupledValue("rhoE"))
{}

Real
InternalEnergyAux::computeValue()
{
    // Compute the velocity:
    Real _vel = _rhou[_qp] / _rho[_qp];
    
    // Return internal energy:
    return _rhoE[_qp] - 0.5*_rho[_qp]*_vel*_vel;
}
