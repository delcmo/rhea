#include "ComputeICsRadHydro.h"

template<>
InputParameters validParams<ComputeICsRadHydro>()
{
  InputParameters params = validParams<UserObject>();

  // Physical properties
  params.addRequiredParam<Real>("speed_of_light", "Speed of light");
  params.addRequiredParam<Real>("a", "Boltzman constant");
  // Non-dimensionalized numbers
  params.addRequiredParam<Real>("P", "Ratio of radiant energy to material energy");
  params.addRequiredParam<Real>("Mach_inlet", "Inlet Mach number");
  // Pre-schock conditions:
  params.addRequiredParam<Real>("rho_hat_0", "Non-dimensionalyzed pre-shock density");
  params.addRequiredParam<Real>("T_hat_0", "Non-dimensionalyzed pre-shock temperature");
  // Equation of state
  params.addRequiredParam<UserObjectName>("eos", "parameters for eos.");

  return params;
}

ComputeICsRadHydro::ComputeICsRadHydro(const std::string & name, InputParameters parameters) :
    GeneralUserObject(name, parameters),
    // Physical properties
    _sp(getParam<Real>("speed_of_light")),
    _a(getParam<Real>("a")),
    // Non-dimensionalized numbers
    _P(getParam<Real>("P")),
    _Mach_inlet(getParam<Real>("Mach_inlet")),
    // Pre-shock parameters
    _rho_hat_pre(getParam<Real>("rho_hat_0")),
    _T_hat_pre(getParam<Real>("T_hat_0")),
    // User Objects
    _eos(getUserObject<IdealGasEquationOfState>("eos"))
{
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
  _T_hat_post = T_1*_T_hat_pre;
  _rho_hat_post = rho_1*_rho_hat_pre;

  /// Compute other pre and post shock parameters ///
  // Compute vel_hat_pre and vel_hat_post
  Real a_hat_0 = std::sqrt(_a*_T_hat_pre*_T_hat_pre*_T_hat_pre*_T_hat_pre/(_rho_hat_pre*_P));
  _vel_hat_pre = _Mach_inlet*a_hat_0;
  _vel_hat_post = _rho_hat_pre*_vel_hat_pre/_rho_hat_post;

  // Compute eps_hat_pre and eps_hat_post
  _eps_hat_pre = _a*_T_hat_pre*_T_hat_pre*_T_hat_pre*_T_hat_pre;
  _eps_hat_post = _a*_T_hat_post*_T_hat_post*_T_hat_post*_T_hat_post;

  // Compute and output the pre- and post-shock density values
  std::cout<<"--------------------------------------------------------------"<<std::endl;
  std::cout<<"Pre-schock density value: "<< _rho_hat_pre << std::endl;
  std::cout<<"Post-schock density value: "<< _rho_hat_post << std::endl;

  // Compute and output the pre- and post-shock momentum values
  std::cout<<"--------------------------------------------------------------"<<std::endl;
  std::cout<<"Pre-schock momentum value: "<< _rho_hat_pre*_vel_hat_pre << std::endl;
  std::cout<<"Post-schock momentum value: "<< _rho_hat_post*_vel_hat_post << std::endl;

  // Compute and output the pre- and post-shock material energy values
  Real press_hat_0 = _eos.p_from_T_rho(_T_hat_pre, _rho_hat_pre);
  Real press_hat_1 = _eos.p_from_T_rho(_T_hat_post, _rho_hat_post);
  Real e_hat_0 = _eos.e_from_p_rho(press_hat_0, _rho_hat_pre);
  Real e_hat_1 = _eos.e_from_p_rho(press_hat_1, _rho_hat_post);
  std::cout<<"--------------------------------------------------------------"<<std::endl;
  std::cout<<"Pre-schock material energy value: "<< _rho_hat_pre*(e_hat_0+0.5*_vel_hat_pre*_vel_hat_pre) << std::endl;
  std::cout<<"Post-schock material energy value: "<< _rho_hat_post*(e_hat_1+0.5*_vel_hat_post*_vel_hat_post) << std::endl;

  // Compute and output the pre- and post-shock radiation energy values
  std::cout<<"--------------------------------------------------------------"<<std::endl;
  std::cout<<"Pre-schock radiation energy value: "<< _eps_hat_pre << std::endl;
  std::cout<<"Post-schock radiation energy value: "<< _eps_hat_post << std::endl;

  // Compute and output the pre- and post-shock pressure values
  std::cout<<"--------------------------------------------------------------"<<std::endl;
  std::cout<<"Pre-schock pressure value: "<< press_hat_0 << std::endl;
  std::cout<<"Post-schock pressure value: "<< press_hat_1 << std::endl;
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
