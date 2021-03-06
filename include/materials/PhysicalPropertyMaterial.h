#ifndef PHYSICALPRPERTYMATERIAL_H
#define PHYSICALPRPERTYMATERIAL_H

#include "Material.h"
#include "MaterialProperty.h"
#include "EquationOfState.h"
#include "InputFileSpecifiedICsRadHydro.h"

//Forward Declarations
class PhysicalPropertyMaterial;

template<>
InputParameters validParams<PhysicalPropertyMaterial>();

class PhysicalPropertyMaterial : public Material
{
public:
  PhysicalPropertyMaterial(const InputParameters & parameters);

protected:
  virtual void computeQpProperties();

private:
  // Cross-section types
  enum CrossSectionType
  {
    constant_cs = 0,
    temp_dpt_cs = 1,
    temp_dpt_opacity = 2,
    pure_absorber = 3
  };
  MooseEnum _cs_type;

  // Boolean for diffusion coefficient and dimensional form
  bool _is_diffusion;
  bool _is_dmsl_form;  

  // Coupled variables:
  const VariableValue & _rho;

  // Coupled aux variables:
  const VariableValue & _press;

  // Material properties: cross-section and diffusion.
  MaterialProperty<Real> & _sigma_a;
  MaterialProperty<Real> & _sigma_t;
  MaterialProperty<Real> & _diffusion;

  // Constant cross-sections
  Real _sigma_hat_a;
  Real _sigma_hat_t;

  // Constant opacity
  Real _opacity_a;
  Real _opacity_t;

  // Vector storing parameters for cross sections
  RealVectorValue _sigma_a0;
  RealVectorValue _sigma_t0;

  // Equation of state
  const EquationOfState & _eos;

  // Userobject computing the ICs
  const InputFileSpecifiedICsRadHydro & _ics;
};

#endif // PHYSICALPRPERTYMATERIAL_H