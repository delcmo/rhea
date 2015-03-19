#ifndef PHYSICALPRPERTYMATERIAL_H
#define PHYSICALPRPERTYMATERIAL_H

#include "Material.h"
#include "MaterialProperty.h"
#include "EquationOfState.h"

//Forward Declarations
class PhysicalPropertyMaterial;

template<>
InputParameters validParams<PhysicalPropertyMaterial>();

class PhysicalPropertyMaterial : public Material
{
public:
  PhysicalPropertyMaterial(const std::string & name, InputParameters parameters);

protected:
  virtual void computeQpProperties();

private:
  // Cross-section types
  enum CrossSectionType
  {
    constant = 0,
    temp_dpt = 1
  };
  MooseEnum _cs_type; 

  // Coupled variables:
  VariableValue & _rho;

  // Coupled aux variables:
  VariableValue & _press;

  // Material properties: cross-section and diffusion.
  MaterialProperty<Real> & _sigma_a;
  MaterialProperty<Real> & _sigma_t;
  MaterialProperty<Real> & _diffusion;

  // Boltzman constant and speed of light:
  Real _c;
  Real _a;

  // Non-dimensionalized parameters
  Real _P;
  Real _K;
  Real _SIGMA_A;
  
  // Pre-shock parameters
  Real _rho_hat_pre;
  Real _T_hat_pre;

  // Constant cross-sections
  Real _sigma_hat_a;
  Real _sigma_hat_t;

  // Vector storing parameters for cross sections
  RealVectorValue _sigma_a0;
  RealVectorValue _sigma_t0;

  // Equation of state
  const EquationOfState & _eos;
};

#endif // PHYSICALPRPERTYMATERIAL_H