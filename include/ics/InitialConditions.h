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

#ifndef INITIALCONDITIONS_H
#define INITIALCONDITIONS_H

// MOOSE Includes
#include "InitialCondition.h"
#include "EquationOfState.h"

// Forward Declarations
class InitialConditions;

template<>
InputParameters validParams<InitialConditions>();

/**
 * InitialConditions just returns a constant value.
 */
class InitialConditions : public InitialCondition
{
public:

  /**
   * Constructor: Same as the rest of the MOOSE Objects
   */
  InitialConditions(const std::string & name,
            InputParameters parameters);

  /**
   * The value of the variable at a point.
   *
   * This must be overriden by derived classes.
   */
  virtual Real value(const Point & p);

private:
    // Initial conditions for left and right values:
    Real _rho_left;
    Real _rho_right;
    Real _v_left;
    Real _v_right;
    Real _t_left;
    Real _t_right;
    Real _eps_left;
    Real _eps_right;
    
    // Position of the membrane:
    Real _membrane;
    
    // Name of the variable:
    std::string _name_var;
    
    // Equation of state:
    const EquationOfState & _eos;

};

#endif //INITIALCONDITIONS_H
