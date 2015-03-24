#include "PhysicalPropertyMaterial.h"

template<>
InputParameters validParams<PhysicalPropertyMaterial>()
{
  InputParameters params = validParams<Material>();

  // Boolean for diffusion
  params.addParam<bool>("is_diffusion", true, "boolean to turn the diffusion on/off");
  // Coupled variable
  params.addRequiredCoupledVar("rho", "density");
  // Coupled aux variable
  params.addRequiredCoupledVar("pressure", "pressure");
  // Type of cross-section
  params.addParam<std::string>("cross_section_name", "constant_cs", "Cross-section type used in the simulation.");
  // Boltzman constant and speed of light:
  params.addParam<Real>("speed_of_light", 299.792, "speed of light");
  params.addParam<Real>("a", 1.372e-2, "Boltzman constant");
  // Non-dimensionalized numbers
  params.addParam<Real>("P", "Ratio of radiant energy to material energy");
  params.addParam<Real>("K", "Radiative diffusivity");
  params.addParam<Real>("SIGMA_A", "Non-dimensionalized absorption cross-section");
  // Pre-shock parameters
  params.addParam<Real>("rho_hat_0", "Pre-shock density value");
  params.addParam<Real>("T_hat_0", "Pre-shock temperature value");
  // Parameters used in the computation of the cross sections
  params.addParam<RealVectorValue>("sigma_t0", "parameters (a,b,c) for total cross-section (a and c are unitless): a / (b + T^c)");
  params.addParam<RealVectorValue>("sigma_a0", "parameters (a,b,c) for absorption cross-section (a and c are unitless): a / (b + T^c)");
  // Userobject:
  params.addRequiredParam<UserObjectName>("eos", "Equation of state");

  return params;
}

PhysicalPropertyMaterial::PhysicalPropertyMaterial(const std::string & name, InputParameters params) :
    Material(name, params),
    // Declare viscosity types
    _cs_type("constant_cs temp_dpt_cs temp_dpt_opacity invalid", getParam<std::string>("cross_section_name")),
    // Boolean for diffusion
    _is_diffusion(getParam<bool>("is_diffusion")),
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
    _P(isParamValid("P") ? getParam<Real>("P") : 1.),
    _K(isParamValid("K") ? getParam<Real>("K") : 1.),
    _SIGMA_A(isParamValid("SIGMA_A") ? getParam<Real>("SIGMA_A") : 0.),
    // Pre-shock parameters
    _rho_hat_pre(isParamValid("rho_hat_0") ? getParam<Real>("rho_hat_0") : 1.),
    _T_hat_pre(isParamValid("T_hat_0") ? getParam<Real>("T_hat_0") : 1.),
    // Parameters used in the computation of the cross sections
    _sigma_a0(isParamValid("sigma_a0") ? getParam<RealVectorValue>("sigma_a0") : (1., 0., 0.)),
    _sigma_t0(isParamValid("sigma_t0") ? getParam<RealVectorValue>("sigma_t0") : (1., 0., 0.)),
    // Equation of state
    _eos(getUserObject<EquationOfState>("eos"))
{
  // Computed the constant total and absorption cross-sections from the initial conditions
  Real a_hat_0 = std::sqrt(_a*_T_hat_pre*_T_hat_pre*_T_hat_pre*_T_hat_pre/(_rho_hat_pre*_P));
  if (_cs_type == constant_cs)
  {
    if ( !params.isParamValid("P") || !params.isParamValid("K") || !params.isParamValid("SIGMA_A") )
      mooseError("'"<<this->name()<<"': the cross section are constant but valid input parameters are not provided: P, K and/or SIGMA_A.");

    if ( !params.isParamValid("rho_hat_0") || !params.isParamValid("T_hat_0") )
      mooseError("'"<<this->name()<<"': the cross section are constant but valid input parameters are not provided: rho_hat_0 and/or T_hat_0.");

    _sigma_hat_t = _c/(3.*_K*a_hat_0);
    _sigma_hat_a = _SIGMA_A*a_hat_0/_c;
  }

  // Check for valid parameters if temperature-dependent cross section
  if (_cs_type == temp_dpt_cs)
  {
    if (!params.isParamValid("sigma_a0") || !params.isParamValid("sigma_t0"))
      mooseError("'"<<this->name()<<"': the cross section are temperature dependent but valid input parameters are not provided: sigma_a0 and/or sigma_t0.");

    _sigma_hat_t = _sigma_t0(0)*a_hat_0/_c;
    _sigma_hat_a = _sigma_a0(0)*a_hat_0/_c;
  }

  // Check for valid parameters if opacity
  if (_cs_type == temp_dpt_opacity)
  {
    if (!params.isParamValid("sigma_a0") || !params.isParamValid("sigma_t0") || !isCoupled("rho"))
      mooseError("'"<<this->name()<<"': the cross sections are computed from the opacity but valid input parameters are not provided: sigma_a0, sigma_t0 and/or density.");

    _opacity_a = _sigma_a0(0)*a_hat_0/_c;
    _opacity_t = _sigma_t0(0)*a_hat_0/_c;
  }
}

void
PhysicalPropertyMaterial::computeQpProperties()
{
  // Cross sections:
  Real temp = 0.;
  switch (_cs_type)
  {
    case constant_cs: // constant cross section
      _sigma_t[_qp] = _sigma_hat_t;
      _sigma_a[_qp] = _sigma_hat_a;
      break;
    case temp_dpt_cs: // temperature-dependent cross section: \sigma_x = \sigma_hat_x / (\sigma_x(1) + T^\sigma_x(2))
      temp = _eos.temperature_from_p_rho(_press[_qp], _rho[_qp]);
      _sigma_t[_qp] = _sigma_hat_t / (_sigma_t0(1) + std::pow(temp, _sigma_t0(2)));
      _sigma_a[_qp] = _sigma_hat_a / (_sigma_a0(1) + std::pow(temp, _sigma_a0(2)));
      break;
    case temp_dpt_opacity: // temperature-dependent cross section: \sigma_x = rho * _opacity_x / (\sigma_x(1) + T^\sigma_x(2))
      temp = _eos.temperature_from_p_rho(_press[_qp], _rho[_qp]);      
      _sigma_t[_qp] = _rho[_qp] * _opacity_t / (_sigma_t0(1) + std::pow(temp, _sigma_t0(2)));
      _sigma_a[_qp] = _rho[_qp] * _opacity_a / (_sigma_a0(1) + std::pow(temp, _sigma_a0(2)));
      break;
    default:
      mooseError("In the function '"<<this->name()<<"', the cross-section type "<< _cs_type <<" is not implemented.");
      break;
  }

  // Diffusion coefficient
  if (_is_diffusion)
    _diffusion[_qp] = _c / (3*_sigma_t[_qp]);
  else
    _diffusion[_qp] = 0.;

  // Check
  if (_sigma_a[_qp]>_sigma_t[_qp])
    mooseError("The absoprtion cross section ("<<_sigma_t[_qp]<< ") coefficient computed in '" << this->name() << "' is larger than the total cross section ("<<_sigma_t[_qp]<<")");
  if (_diffusion[_qp]<0)
    mooseError("The diffusion coefficient computed in '" << this->name() << "' is negative.");
}