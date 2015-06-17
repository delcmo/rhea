#include "PhysicalPropertyMaterial.h"

template<>
InputParameters validParams<PhysicalPropertyMaterial>()
{
  InputParameters params = validParams<Material>();

  // Booleans for diffusion and dimensional form
  params.addParam<bool>("is_diffusion", true, "boolean to turn the diffusion on/off");
  params.addParam<bool>("is_dimensional_form", true, "boolean to solve the momentum equation in a dimensional form");
  // Coupled variable
  params.addRequiredCoupledVar("rho", "density");
  // Coupled aux variable
  params.addRequiredCoupledVar("pressure", "pressure");
  // Type of cross-section
  params.addParam<std::string>("cross_section_name", "constant_cs", "Cross-section type used in the simulation.");
  // Parameters used in the computation of the cross sections
  params.addParam<RealVectorValue>("sigma_t0", "parameters (a,b,c) for total cross-section (a and c are unitless): a / (b + T^c)");
  params.addParam<RealVectorValue>("sigma_a0", "parameters (a,b,c) for absorption cross-section (a and c are unitless): a / (b + T^c)");
  // Userobject:
  params.addRequiredParam<UserObjectName>("eos", "Equation of state");
  // Userobject computing the ICs
  params.addRequiredParam<UserObjectName>("ics", "parameters for ics.");

  return params;
}

PhysicalPropertyMaterial::PhysicalPropertyMaterial(const std::string & name, InputParameters params) :
    Material(name, params),
    // Declare viscosity types
    _cs_type("constant_cs temp_dpt_cs temp_dpt_opacity pure_absorber invalid", getParam<std::string>("cross_section_name")),
    // Booleans
    _is_diffusion(getParam<bool>("is_diffusion")),
    _is_dmsl_form(getParam<bool>("is_dimensional_form")),
    // Coupled variables:
    _rho(coupledValue("rho")),
    // Coupled aux variables
    _press(coupledValue("pressure")),
    // Declare material properties: cross-section and diffusion.
    _sigma_a(declareProperty<Real>("sigma_a")),
    _sigma_t(declareProperty<Real>("sigma_t")),
    _diffusion(declareProperty<Real>("diffusion")),
    // Parameters used in the computation of the cross sections
    _sigma_a0(isParamValid("sigma_a0") ? getParam<RealVectorValue>("sigma_a0") : (1., 0., 0.)),
    _sigma_t0(isParamValid("sigma_t0") ? getParam<RealVectorValue>("sigma_t0") : (1., 0., 0.)),
    // Equation of state
    _eos(getUserObject<EquationOfState>("eos")),
    // Userobject computing the ICs
    _ics(getUserObject<ComputeICsRadHydro>("ics"))
{
  // Computed the constant total and absorption cross-sections from the initial conditions
  Real a_hat_0 = std::sqrt(_ics.a()*_ics.T_hat_pre()*_ics.T_hat_pre()*_ics.T_hat_pre()*_ics.T_hat_pre()/(_ics.rho_hat_pre()*_ics.P()));
  if (_cs_type == constant_cs)
  {
    _sigma_hat_t = _ics.c()/(3.*_ics.K()*a_hat_0);
    _sigma_hat_a = _ics.SIGMA_A()*a_hat_0/_ics.c();

    // Output the constant total and absorption cross sections
    std::cout<<"Constant total cross section: "<< _sigma_hat_t << std::endl;
    std::cout<<"Constant absorption cross section: "<< _sigma_hat_a << std::endl;
    std::cout<<"--------------------------------------------------------------"<<std::endl;    
  }

  // Temperature-dependent cross section
  if (_cs_type == temp_dpt_cs)
  {
    if (!params.isParamValid("sigma_a0") || !params.isParamValid("sigma_t0"))
      mooseError("'"<<this->name()<<"': the cross section are temperature dependent but valid input parameters are not provided: sigma_a0 and/or sigma_t0.");

    _sigma_hat_t = _sigma_t0(0)*a_hat_0/_ics.c();
    _sigma_hat_a = _sigma_a0(0)*a_hat_0/_ics.c();
  }

  // Constant opacity
  if (_cs_type == temp_dpt_opacity)
  {
    if (!params.isParamValid("sigma_a0") || !params.isParamValid("sigma_t0") || !isCoupled("rho"))
      mooseError("'"<<this->name()<<"': the cross sections are computed from the opacity but valid input parameters are not provided: sigma_a0, sigma_t0 and/or density.");

    _opacity_a = _sigma_a0(0)*std::pow(_ics.T_hat_pre(), 3.5)/_ics.rho_hat_pre();
    _opacity_t = _sigma_t0(0)*std::pow(_ics.T_hat_pre(), 3.5)/_ics.rho_hat_pre();
  }

  // Constant cross section for pure aborber material
  if (_cs_type == pure_absorber)
  {
    Real C0 = std::sqrt(_ics.SIGMA_A()*3.*_ics.K());
    _sigma_hat_a = _ics.SIGMA_A()/C0;
    _sigma_hat_t = _sigma_hat_a;
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
      _sigma_t[_qp] = _is_dmsl_form ? _sigma_hat_t : 1.;
      _sigma_a[_qp] = _is_dmsl_form ? _sigma_hat_a : 1.;
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
    case pure_absorber:
      _sigma_t[_qp] = _is_dmsl_form ? _sigma_hat_t : 1.;
      _sigma_a[_qp] = _is_dmsl_form ? _sigma_hat_a : 1.;
      break;
    default:
      mooseError("In the function '"<<this->name()<<"', the cross-section type "<< _cs_type <<" is not implemented.");
      break;
  }

  // Diffusion coefficient
  if (_is_diffusion)
    _diffusion[_qp] = _is_dmsl_form ? _ics.c() / (3*_sigma_t[_qp]) : 1. / _sigma_t[_qp];
  else
    _diffusion[_qp] = 0.;

  // Check
  if (_sigma_a[_qp]>_sigma_t[_qp])
    mooseError("The absoprtion cross section ("<<_sigma_a[_qp]<< ") coefficient computed in '" << this->name() << "' is larger than the total cross section ("<<_sigma_t[_qp]<<")");
  if (_diffusion[_qp]<0)
    mooseError("The diffusion coefficient computed in '" << this->name() << "' is negative.");
}