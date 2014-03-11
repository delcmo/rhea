#ifndef COMPUTEMATERIALS_H
#define COMPUTEMATERIALS_H

#include "Material.h"
#include "MaterialProperty.h"
#include "EquationOfState.h"

//Forward Declarations
class ComputeMaterials;

template<>
InputParameters validParams<ComputeMaterials>();

class ComputeMaterials : public Material
{
public:
  ComputeMaterials(const std::string & name, InputParameters parameters);

protected:
  virtual void computeQpProperties();

private:
    // Cross-section types
    enum CrossSectionType
    {
        CONSTANT = 0,
        TYPE1 = 1
    };
    
    // Viscosity types
    enum ViscosityType
    {
        FIRST_ORDER = 0,
        ENTROPY = 1
    };
    
    // Cross-section:
    std::string _cs_name;
    MooseEnum _cs_type;
    
    // Artificial diffusion name
    std::string _visc_name;
    
    // Aritifical diffusion type
    MooseEnum _visc_type;
    
    // Variables:
    VariableValue & _vel;
    VariableValue & _vel_old;
    VariableValue & _rho;
    VariableValue & _rho_old;
    VariableGradient & _grad_rho;
    VariableGradient & _grad_rho_old;
    VariableValue & _pressure;
    VariableValue & _pressure_old;
    VariableGradient & _grad_press;
    VariableGradient & _grad_press_old;
    VariableValue & _epsilon;
    VariableValue & _epsilon_old;
    VariableGradient & _grad_eps;
    VariableGradient & _grad_eps_old;
    
    // Jumps:
    VariableValue & _jump_press;
    VariableValue & _jump_dens;
    
    // Material properties: cross-section and diffusion.
    MaterialProperty<Real> & _sigma_a;
    MaterialProperty<Real> & _sigma_t;
    MaterialProperty<Real> & _diffusion;
    
    // Material properties: viscosity coefficients.
    MaterialProperty<Real> & _mu;
    MaterialProperty<Real> & _mu_max;
    MaterialProperty<Real> & _kappa;
    MaterialProperty<Real> & _kappa_max;
    
    // Material constants:
    Real _c;
    Real _a;
    Real _sigma_a0;
    Real _sigma_t0;

    // Multiplicative coefficient for viscosity:
    double _Ce;
    
    // UserObject: equation of state
    const EquationOfState & _eos;
    
    // Name of the posprocessors for pressure and velocity:
    std::string _epsilon_pps_name;
    std::string _velocity_pps_name;
};

#endif //COMPUTEMATERIALS_H
