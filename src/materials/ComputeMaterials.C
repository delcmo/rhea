#include "ComputeMaterials.h"

template<>
InputParameters validParams<ComputeMaterials>()
{
  InputParameters params = validParams<Material>();
    // Type of viscosity:
    params.addParam<std::string>("viscosity_name", "FIRST_ORDER", "Name of the viscosity definition to use: set to FIRST ORDER by default.");
    // Boolean for cross-section:
    params.addParam<bool>("TempDepedentCrossSection", false, "boolean to make the cross-sections temperature dependent. False by default.");
    // Coupled variables:
    params.addRequiredCoupledVar("velocity", "velocity");
    params.addRequiredCoupledVar("density", "density");
    params.addRequiredCoupledVar("pressure", "pressure");
    params.addCoupledVar("epsilon", "epsilon");
    params.addCoupledVar("jump", "jumps of the pressure gradient");
    // Material constants:
    params.addRequiredParam<Real>("speed_of_light", "speed of light");
    params.addRequiredParam<Real>("a", "a");
    params.addRequiredParam<Real>("sigma_a0", "absorption cross-section");
    params.addRequiredParam<Real>("sigma_t0", "total cross-section");
    // Cconstant parameter:
    params.addParam<double>("Ce", 1., "Coefficient for viscosity");
    // Userobject:
    params.addRequiredParam<UserObjectName>("eos", "Equation of state");
    // PPS names:
    params.addParam<std::string>("epsilon_PPS_name", "none", "name of the pps for radiation");
    params.addParam<std::string>("velocity_PPS_name", "none", "name of the pps for velocity");
    return params;
}

ComputeMaterials::ComputeMaterials(const std::string & name, InputParameters parameters) :
    Material(name, parameters),
    // Boolean for temperature dependent cross-section:
    _tempDependentCS(getParam<bool>("TempDepedentCrossSection")),
    // Declare viscosity types:
    _visc_name(getParam<std::string>("viscosity_name")),
    _visc_type("FIRST_ORDER, ENTROPY, INVALID", "INVALID"),
    // Variables:
    _vel(coupledValue("velocity")),
    _vel_old(coupledValueOld("velocity")),
    _rho(coupledValue("density")),
    _rho_old(coupledValueOld("density")),
    _grad_rho(coupledGradient("density")),
    _grad_rho_old(coupledGradientOld("density")),
    _pressure(coupledValue("pressure")),
    _pressure_old(coupledValueOld("pressure")),
    _grad_press(coupledGradient("pressure")),
    _grad_press_old(coupledGradientOld("pressure")),
    _epsilon(isCoupled("epsilon") ?  coupledValue("epsilon") : _zero),
    _epsilon_old(isCoupled("epsilon") ? coupledValueOld("epsilon") : _zero),
    _grad_eps(isCoupled("epsilon") ? coupledGradient("epsilon") : _grad_zero),
    _grad_eps_old(isCoupled("epsilon") ? coupledGradientOld("epsilon") : _grad_zero),
    // Jumps:
    _jump(isCoupled("jump") ? coupledValue("jump") : _zero),
    // Declare material properties: cross-section and diffusion.
    _sigma_a(declareProperty<Real>("sigma_a")),
    _sigma_t(declareProperty<Real>("sigma_t")),
    _diffusion(declareProperty<Real>("diffusion")),
    // Declare material properties for viscosity coefficients.
    _mu(declareProperty<Real>("mu")),
    _mu_max(declareProperty<Real>("mu_max")),
    _kappa(declareProperty<Real>("kappa")),
    _kappa_max(declareProperty<Real>("kappa_max")),
    // Material constants:
    _c(getParam<Real>("speed_of_light")),
    _a(getParam<Real>("a")),
    _sigma_a0(getParam<Real>("sigma_a0")),
    _sigma_t0(getParam<Real>("sigma_t0")),
    // Get parameter Ce
    _Ce(getParam<double>("Ce")),
    // UserObject:
    _eos(getUserObject<EquationOfState>("eos")),
    // PPS name:
    _epsilon_pps_name(getParam<std::string>("epsilon_PPS_name")),
    _velocity_pps_name(getParam<std::string>("velocity_PPS_name"))
{
    _visc_type = _visc_name;
    if (_Ce < 0.)
        mooseError("The coefficient Ce has to be positive.");
}

void
ComputeMaterials::computeQpProperties()
{
    // Material cross-sections and diffusion:
    _sigma_a[_qp] = _sigma_a0;
    _sigma_t[_qp] = _sigma_t0;
    _diffusion[_qp] = _c / (3*_sigma_t[_qp]);
    
    if (_tempDependentCS) {
        Real _temp = _eos.temperature_from_p_rho(_pressure[_qp], _rho[_qp]);
        _diffusion[_qp] = _c / (3*_sigma_t[_qp]);
        _sigma_a[_qp] = _sigma_a0;
        _sigma_t[_qp] = _sigma_t0;
    }
    
    // Determine h (length used in definition of first and second order viscosities):
    Real _h = _current_elem->hmin();
    
    // Compute first order viscosity:
    //std::cout<<"press="<<_pressure[_qp]<<std::endl;
    //std::cout<<"rho="<<_rho[_qp]<<std::endl;
    Real _c = std::sqrt(_eos.c2_from_p_rho(_rho[_qp], _pressure[_qp]));
    //std::cout<<"c="<<_c<<std::endl;
    _mu_max[_qp] = 0.5*_h*(_vel[_qp] + _c);
    _kappa_max[_qp] = 0.5*_h*(_vel[_qp] + _c);
    
    // Get postprocessor values and compute the normalization factor:
    Real _eps = std::sqrt(std::numeric_limits<Real>::min());
    Real _epsilon_pps = std::max(getPostprocessorValueByName(_epsilon_pps_name), _eps);
    Real _vel_pps = std::max(getPostprocessorValueByName(_velocity_pps_name), _eps);
    Real _norm_mu = 0.5*_rho[_qp]*_c*_c;
    Real _norm_kappa = 0.5*std::min(_rho[_qp]*_c*_c, _rho[_qp]*_c*_vel_pps);
    
    // Switch statement over viscosity type:
    Real _kappa_e = 0.; Real _mu_e = 0.; Real _jump_value = 0.; Real _residual = 0;
    Real _DP = 0.; Real _Drho = 0.; Real _Deps = 0.; Real _Dstt = 0.;
    Real _vel_half = 0.5*(_vel[_qp]+_vel_old[_qp]);
    switch (_visc_type) {
        case FIRST_ORDER:
            _mu[_qp] = _mu_max[_qp];
            _kappa[_qp] = _kappa_max[_qp];
            break;
        case ENTROPY:
            // Compute the jump:
            _jump_value = (std::fabs(_vel[_qp])+_c)*_jump[_qp];
            
            // Compute the pressure residual:
            _Dstt = 0.5*_vel_half*(_grad_press[_qp](0)+_grad_press_old[_qp](0));
            _DP = (_pressure[_qp]-_pressure_old[_qp])/_dt + _Dstt;
            
            // Compute the density residual:
            _Dstt = 0.5*_vel_half*(_grad_rho[_qp](0)+_grad_rho_old[_qp](0));
            _Drho = (_rho[_qp]-_rho_old[_qp])/_dt + _Dstt;
            
            // Compute radiation residual:
            _Dstt = 0.5*_vel_half*(_grad_eps[_qp](0)+_grad_eps_old[_qp](0));
            _Deps = 0;//(_epsilon[_qp]-_epsilon_old[_qp])/_dt + _Dstt;
            
            // Compute mu_e and kappa_e:
            _residual = std::fabs( _DP + _c*_c*_Drho + _Deps );
            _mu_e = _Ce*_h*_h*std::max(_residual, 2*_jump_value) / _norm_mu;
            _kappa_e = _Ce*_h*_h*std::max(_residual, 2*_jump_value) / _norm_kappa;
            
            // Get the value of the viscosity coefficients:
            _kappa[_qp] = std::min( _kappa_max[_qp], _kappa_e );
            _mu[_qp] = std::min( _mu_max[_qp], _mu_e );
            break;
        default:
            mooseError("The viscosity type entered in the input file is not implemented.");
            break;
    }
}
