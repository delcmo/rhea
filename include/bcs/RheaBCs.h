#ifndef RHEABCS_H
#define RHEABCS_H

#include "IntegratedBC.h"
#include "EquationOfState.h"

// Forward Declarations
class RheaBCs;
class EquationOfState;

template<>
InputParameters validParams<RheaBCs>();

class RheaBCs : public IntegratedBC
{

public:
  RheaBCs(const std::string & name, InputParameters parameters);

  virtual ~RheaBCs(){}

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned jvar);

    enum EFlowEquationType
    {
        CONTINUITY = 0,
        MOMENTUM = 1,
        ENERGY = 2,
        RADIATION = 3
    };
    
    // Boolean for coupled velocity
    bool _isVel;

    // Eqn. name to be read from input file
    std::string _eqn_name;
    
    // which equation (mass/momentum/energy) this BC is acting on
    MooseEnum _eqn_type;
    
    // Coupled aux variables
    VariableValue & _vel;
    VariableValue & _rho;
    VariableValue & _pressure;
    VariableValue & _epsilon;
    VariableValue & _vel_old;
    VariableValue & _rho_old;
    VariableValue & _pressure_old;
    
    // Specified pressure
    Real _p_bc;
    
    // Specified temperature
    Real _T_bc;
    
    // Specified velocity:
    Real _v_bc;
    
    // Specified radiation:
    Real _eps_left;
    Real _eps_right;
    Real _c_light;
    
    // Material property:
    MaterialProperty<Real> & _D;
    
    // Equation of state
    const EquationOfState & _eos;
};

#endif // RHEABCS_H

