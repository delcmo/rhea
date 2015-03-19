#include "PhysicalPropertyMaterial.h"

template<>
InputParameters validParams<PhysicalPropertyMaterial>()
{
  InputParameters params = validParams<Material>();

  // Coupled variable
  params.addRequiredCoupledVar("rho", "density");
  // Coupled aux variable
  params.addRequiredCoupledVar("pressure", "pressure");
  // Type of cross-section
  params.addParam<std::string>("cross_section_name", "constant", "Cross-section type used in the simulation.");
  // Boltzman constant and speed of light:
  params.addParam<Real>("speed_of_light", 299.792, "speed of light");
  params.addParam<Real>("a", 1.372e-2, "Boltzman constant");
  // Non-dimensionalized numbers
  params.addParam<Real>("P", -1., "Ratio of radiant energy to material energy");
  params.addParam<Real>("K", -1., "Radiative diffusivity");
  params.addParam<Real>("SIGMA_A", -1., "Non-dimensionalized absorption cross-section");
  // Pre-shock parameters
  params.addParam<Real>("rho_hat_0", -1., "Pre-shock density value");
  params.addParam<Real>("T_hat_0", -1., "Pre-shock temperature value");
  // Parameters used in the computation of the cross sections  
  params.addParam<RealVectorValue>("sigma_a0", (-1., 0., 0), "absorption cross-section coefficient (sigma_a0, sigma_a1, t_a)");
  params.addParam<RealVectorValue>("sigma_t0", (-1., 0., 0), "total cross-section coefficient (sigma_t0, sigma_t1, n_t)");
  // Userobject:
  params.addRequiredParam<UserObjectName>("eos", "Equation of state");

  return params;
}

PhysicalPropertyMaterial::PhysicalPropertyMaterial(const std::string & name, InputParameters parameters) :
    Material(name, parameters),
    // Declare viscosity types
    _cs_type("constant temp_dpt invalid", getParam<std::string>("cross_section_name")),
    // Coupled variables:
    _rho(coupledValue("rho")),
    // Coupled aux variables
    _press(coupledValue("pressure")),
    // Declare material properties: cross-section and diffusion.
    _sigma_a(declareProperty<Real>("sigma_a")),
    _sigma_t(declareProperty<Real>("sigma_t")),
    _diffusion(declareProperty<Real>("diffusion")),
    // Boltzman constant and speed of light:
    _c(getParam<Real>("speed_of_light")),
    _a(getParam<Real>("a")),
    // Non-dimensionalized numbers
    _P(getParam<Real>("P")),
    _K(getParam<Real>("K")),
    _SIGMA_A(getParam<Real>("SIGMA_A")),
    // Pre-shock parameters
    _rho_hat_pre(getParam<Real>("rho_hat_0")),
    _T_hat_pre(getParam<Real>("T_hat_0")),
    // Parameters used in the computation of the cross sections
    _sigma_a0(getParam<RealVectorValue>("sigma_a0")),
    _sigma_t0(getParam<RealVectorValue>("sigma_t0")),
    // Equation of state
    _eos(getUserObject<EquationOfState>("eos"))
{
  // Check consistency of the physical parameters
  if (_sigma_a0(0)>_sigma_t0(0))
    mooseError("The value of the absorption cross section provided in '" << this->name() << "' is larger than the total cross-section value.");

  // Computed the constant total and absorption cross-sections from the initial conditions
  if (_cs_type == constant)
  {
    if (_P<0 || _K<0 || _SIGMA_A<0 || _rho_hat_pre<0 || _T_hat_pre<0)
      mooseError("'"<<this->name()<<"': the cross section are constant but valid input parameters are not provided.");

    Real a_hat_0 = std::sqrt(_a*_T_hat_pre*_T_hat_pre*_T_hat_pre*_T_hat_pre/(_rho_hat_pre*_P));
    _sigma_hat_t = _c/(3.*_K*a_hat_0);
    _sigma_hat_a = _SIGMA_A*a_hat_0/_c;
  }

  // Check for valid parameters if temperature-dependent cross section
  if (_cs_type == temp_dpt)
    if (_sigma_a0(0)<0 || _sigma_t0(0)<0)
      mooseError("'"<<this->name()<<"': the cross section are temperature dependent but valid input parameters are not provided.");
}

void
PhysicalPropertyMaterial::computeQpProperties()
{
  // Cross sections:
  Real temp = 0.;
  switch (_cs_type)
  {
    case constant: // constant cross section
      _sigma_t[_qp] = _sigma_hat_t;
      _sigma_a[_qp] = _sigma_hat_a;
      break;
    case temp_dpt: // temperature-dependent cross section: \sigma_x = \sigma_x(0) / (\sigma_x(1) + T^\sigma_x(2)) 
      temp = _eos.temperature_from_p_rho(_press[_qp], _rho[_qp]);
      _sigma_t[_qp] = _sigma_t0(0) / (_sigma_t0(1) + std::pow(temp, _sigma_t0(2)));
      _sigma_a[_qp] = _sigma_a0(0) / (_sigma_a0(1) + std::pow(temp, _sigma_a0(2)));
      break;
    default:
      mooseError("'"<<this->name()<<"':The cross-section type is not implemented.");
      break;
  }

  // Diffusion coefficient
  _diffusion[_qp] = _c / (3*_sigma_t[_qp]);

  // Check
  if (_sigma_a[_qp]>_sigma_t[_qp])
    mooseError("The absoprtion cross section ("<<_sigma_t[_qp]<< ") coefficient computed in '" << this->name() << "' is larger than the total cross section ("<<_sigma_t[_qp]<<")");
  if (_diffusion[_qp]<0)
    mooseError("The diffusion coefficient computed in '" << this->name() << "' is negative.");
}