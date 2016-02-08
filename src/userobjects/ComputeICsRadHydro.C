#include "ComputeICsRadHydro.h"

/** The algorithm used to compute the pre and post shock values si taken from Jarrod Edward's dissertation on page 82 **/

template<>
InputParameters validParams<ComputeICsRadHydro>()
{
  InputParameters params = validParams<UserObject>();

  // Boolean
  params.addParam<bool>("is_dimensional_form", true, "boolean to solve the momentum equation in a dimensional form");
  // Physical properties
  params.addRequiredParam<Real>("speed_of_light", "Speed of light");
  params.addRequiredParam<Real>("a", "Boltzman constant");
  // Non-dimensionalized numbers
  params.addRequiredParam<Real>("P", "Ratio of radiant energy to material energy");
  params.addParam<Real>("K", "Radiative diffusivity");
  params.addParam<Real>("SIGMA_A", "Non-dimensionalized absorption cross-section");
  params.addParam<Real>("C", "Non-dimensionalized ratio of the speed of light to the material speed of sound");
  params.addRequiredParam<Real>("Mach_inlet", "Inlet Mach number");
  // Equation of state
  params.addRequiredParam<UserObjectName>("eos", "parameters for eos.");
  // Execute on set to initial by default
  params.set<MultiMooseEnum>("execute_on") = "initial";

  return params;
}

ComputeICsRadHydro::ComputeICsRadHydro(const InputParameters & parameters) :
    GeneralUserObject(parameters),
    // Boolean
    _is_dmsl_form(getParam<bool>("is_dimensional_form")),
    // Physical properties
    _sp(getParam<Real>("speed_of_light")),
    _a(getParam<Real>("a")),
    // Non-dimensionalized numbers
    _P(getParam<Real>("P")),
    _Mach_inlet(getParam<Real>("Mach_inlet")),
    _K(isParamValid("K") ? getParam<Real>("K") : 1.),
    _SIGMA_A(isParamValid("SIGMA_A") ? getParam<Real>("SIGMA_A") : 0.),
    _C(isParamValid("C") ? getParam<Real>("C") : 1.),
    // User Objects
    _eos(getUserObject<IdealGasEquationOfState>("eos"))
{
  /// Compute rho_0 and T_0
  Real a_hat_0 = _sp/_C;
  Real T_hat_0 = a_hat_0*a_hat_0/(_eos.Cv()*_eos.gamma()*(_eos.gamma()-1));
  Real rho_hat_0 = _a*T_hat_0*T_hat_0*T_hat_0*T_hat_0/(a_hat_0*a_hat_0*_P);

/// Solve for the post-chock temperature T_hat_post and compute the corresponding post-shock density rho_hat_post ///
  // Initialyze value of T_post
  Real T_1 = 1.-_eos.gamma()+2.*_eos.gamma()*_Mach_inlet*_Mach_inlet;
  T_1 *= (2.+(_eos.gamma()-1)*_Mach_inlet*_Mach_inlet)*(2.+(_eos.gamma()-1)*_Mach_inlet*_Mach_inlet);
  T_1 *= 1./((_eos.gamma()+1)*(_eos.gamma()+1)*_Mach_inlet*_Mach_inlet);
  if (_P > 1)
    T_1 = std::pow(8./7.*(_Mach_inlet*_Mach_inlet*9./(4.*_P)-1.) ,0.25);

  // Solve for T_post (Newton solve)
  Real delta_T_1 = 1.;
  while (std::fabs(delta_T_1)>1.e-6)
  {
    Real T_14 = T_1*T_1*T_1*T_1;
    Real f_1 = 3.*(_eos.gamma()+1)*(T_1-1.)-_P*_eos.gamma()*(_eos.gamma()-1.)*(7.+T_14);
    Real f_2 = 12.*(_eos.gamma()-1.)*(_eos.gamma()-1.)*T_1*(3.+_eos.gamma()*_P*(1.+7.*T_14));
    Real rho_1 = (f_1+std::sqrt(f_1*f_1+f_2)) / (6.*(_eos.gamma()-1.)*T_1);
    Real res = 3.*rho_1*(rho_1*T_1-1.);
    res += _eos.gamma()*_P*rho_1*(T_14-1.);
    res -= 3.*_eos.gamma()*(rho_1-1)*_Mach_inlet*_Mach_inlet;
    Real res_prime = 3.*rho_1*rho_1+4.*_eos.gamma()*_P*rho_1*T_1*T_1*T_1;
    delta_T_1 = -res/res_prime;
    T_1 += delta_T_1;

    if (T_1<0)
    {
      T_1=std::fabs(T_1);
      mooseWarning("'"<<this->name()<<"': the post-shock temperature value computed is negative");
    }
  }

  // Compute rho_post
  Real T_14 = T_1*T_1*T_1*T_1;
  Real f_1 = 3.*(_eos.gamma()+1)*(T_1-1.)-_P*_eos.gamma()*(_eos.gamma()-1.)*(7.+T_14);
  Real f_2 = 12.*(_eos.gamma()-1.)*(_eos.gamma()-1.)*T_1*(3.+_eos.gamma()*_P*(1.+7.*T_14));
  Real rho_1 = (f_1+std::sqrt(f_1*f_1+f_2)) / (6.*(_eos.gamma()-1.)*T_1);

  // Compute T_hat_post and rho_hat_post
  _T_hat_pre = _is_dmsl_form ? T_hat_0 : 1.;
  _T_hat_post = _is_dmsl_form ? T_1*T_hat_0 : T_1;
  _rho_hat_pre = _is_dmsl_form ? rho_hat_0 : 1.;
  _rho_hat_post = _is_dmsl_form ? rho_1*rho_hat_0 : rho_1;

  /// Compute other pre and post shock parameters ///
  // Compute vel_hat_pre and vel_hat_post
  _vel_hat_pre = _is_dmsl_form ? _Mach_inlet*a_hat_0 : _Mach_inlet;
  _vel_hat_post = _rho_hat_pre*_vel_hat_pre/_rho_hat_post;

  // Compute eps_hat_pre and eps_hat_post
  _eps_hat_pre = _is_dmsl_form ? _a*_T_hat_pre*_T_hat_pre*_T_hat_pre*_T_hat_pre : _T_hat_pre*_T_hat_pre*_T_hat_pre*_T_hat_pre;
  _eps_hat_post = _is_dmsl_form ? _a*_T_hat_post*_T_hat_post*_T_hat_post*_T_hat_post : _T_hat_post*_T_hat_post*_T_hat_post*_T_hat_post;

  // Compute and output the pre- and post-shock density values
  std::cout<<"--------------------------------------------------------------"<<std::endl;
  std::cout.precision(10);
  std::cout<<"Pre-schock dimensional density value: "<< _rho_hat_pre << std::endl;
  std::cout<<"Post-schock dimensional density value: "<< _rho_hat_post << std::endl;
  std::cout<<"--------------------------------------------------------------"<<std::endl;

  // Compute and output the pre- and post-shock momentum values
  std::cout<<"--------------------------------------------------------------"<<std::endl;
  std::cout<<"Pre-schock dimensional momentum value: "<< _rho_hat_pre*_vel_hat_pre << std::endl;
  std::cout<<"Post-schock dimensional momentum value: "<< _rho_hat_post*_vel_hat_post << std::endl;
  std::cout<<"--------------------------------------------------------------"<<std::endl;

  // Compute and output the pre- and post-shock material energy values
  Real press_hat_pre = _eos.p_from_T_rho(_T_hat_pre, _rho_hat_pre);
  Real press_hat_post = _eos.p_from_T_rho(_T_hat_post, _rho_hat_post);
  Real e_hat_pre = _eos.e_from_p_rho(press_hat_pre, _rho_hat_pre);
  Real e_hat_post = _eos.e_from_p_rho(press_hat_post, _rho_hat_post);
  std::cout<<"--------------------------------------------------------------"<<std::endl;
  std::cout<<"Pre-schock dimensional material energy value: "<< _rho_hat_pre*(e_hat_pre+0.5*_vel_hat_pre*_vel_hat_pre) << std::endl;
  std::cout<<"Post-schock dimensional material energy value: "<< _rho_hat_post*(e_hat_post+0.5*_vel_hat_post*_vel_hat_post) << std::endl;
  std::cout<<"--------------------------------------------------------------"<<std::endl;

  // Compute and output the pre- and post-shock radiation energy values
  std::cout<<"--------------------------------------------------------------"<<std::endl;
  std::cout<<"Pre-schock dimensional radiation energy value: "<< _eps_hat_pre << std::endl;
  std::cout<<"Post-schock dimensional radiation energy value: "<< _eps_hat_post << std::endl;
  std::cout<<"--------------------------------------------------------------"<<std::endl;

  // Compute and output the pre- and post-shock pressure values
  std::cout<<"--------------------------------------------------------------"<<std::endl;
  std::cout<<"Pre-schock dimensional material pressure value: "<< press_hat_pre << std::endl;
  std::cout<<"Post-schock dimensional material pressure value: "<< press_hat_post << std::endl;
  std::cout<<"--------------------------------------------------------------"<<std::endl;

  // Compute and output the pre- and post-shock temperature values
  std::cout<<"--------------------------------------------------------------"<<std::endl;
  std::cout<<"Pre-schock dimensional material temperature value: "<< _T_hat_pre << std::endl;
  std::cout<<"Post-schock dimensional materialtemperature value: "<< _T_hat_post << std::endl;
  std::cout<<"--------------------------------------------------------------"<<std::endl;

  // Compute and output the pre- and post-shock material velocity values
  std::cout<<"--------------------------------------------------------------"<<std::endl;
  std::cout<<"Pre-schock dimensional velocity value: "<< _vel_hat_pre << std::endl;
  std::cout<<"Post-schock dimensional velocity value: "<< _vel_hat_post << std::endl;
  std::cout<<"--------------------------------------------------------------"<<std::endl;

  // Output dimensionless density values
  std::cout<<"--------------------------------------------------------------"<<std::endl;
  std::cout<<"Pre-schock dimensionless density value: "<< 1. << std::endl;
  std::cout<<"Post-schock dimensionless density value: "<< _rho_hat_post/rho_hat_0 << std::endl;
  std::cout<<"--------------------------------------------------------------"<<std::endl;

  // Output dimensionless velocity values
  std::cout<<"--------------------------------------------------------------"<<std::endl;
  std::cout<<"Pre-schock dimensionless velocity value: "<< _vel_hat_pre/a_hat_0 << std::endl;
  std::cout<<"Post-schock dimensionless velocity value: "<< _vel_hat_post/a_hat_0 << std::endl;
  std::cout<<"--------------------------------------------------------------"<<std::endl;

  // Output dimensionless material temperature values
  std::cout<<"--------------------------------------------------------------"<<std::endl;
  std::cout<<"Pre-schock dimensionless material temperature value: "<< 1. << std::endl;
  std::cout<<"Post-schock dimensionless material temperature value: "<< T_1 << std::endl;
  std::cout<<"--------------------------------------------------------------"<<std::endl;

  // Compute and output the pre- and post-shock RED values
  Real red_post = T_1*T_1*T_1*T_1;
  std::cout<<"--------------------------------------------------------------"<<std::endl;
  std::cout<<"Pre-schock dimensionless RED value: "<< 1. << std::endl;
  std::cout<<"Post-schock dimensionless RED value: "<< red_post << std::endl;
  std::cout<<"--------------------------------------------------------------"<<std::endl;

  // Compute and output the pre- and post-shock mach values
  Real mach_hat_pre = _vel_hat_pre / std::sqrt(_eos.c2_from_p_rho(_rho_hat_pre, press_hat_pre, 0.));
  Real mach_hat_post = _vel_hat_post / std::sqrt(_eos.c2_from_p_rho(_rho_hat_post, press_hat_post, 0.));
  std::cout<<"--------------------------------------------------------------"<<std::endl;
  std::cout<<"Pre-schock mach value: "<< mach_hat_pre << std::endl;
  std::cout<<"Post-schock mach value: "<< mach_hat_post << std::endl;
  std::cout<<"--------------------------------------------------------------"<<std::endl;
}

ComputeICsRadHydro::~ComputeICsRadHydro()
{
  // Destructor, empty
}

void
ComputeICsRadHydro::destroy()
{
}