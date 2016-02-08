#include "EntropyViscosityCoefficient.h"
#include "MooseMesh.h"

template<>
InputParameters validParams<EntropyViscosityCoefficient>()
{
  InputParameters params = validParams<Material>();

  // Boolean
  params.addParam<bool>("is_dimensional_form", true, "boolean to solve the equations in a dimensional form");
  // Coupled variables:
  params.addRequiredCoupledVar("rho", "density");
  params.addRequiredCoupledVar("rhou", "momentum");
  params.addCoupledVar("epsilon", "epsilon");
  // Coupled aux variables
  params.addRequiredCoupledVar("pressure", "pressure");
  // Jump values
  params.addCoupledVar("jump_press", "jumps of the pressure gradient");
  params.addCoupledVar("jump_dens", "jumps of the density gradient");
  // Cconstant parameter:
  params.addParam<double>("Cjump", 1., "Coefficient for viscosity");
  params.addParam<bool>("is_first_order_viscosity", false, "Boolean for first order viscosity");
  params.addParam<bool>("use_jumps", true, "Use jumps");
  // Userobject:
  params.addRequiredParam<UserObjectName>("eos", "Equation of state");
  // Userobject computing the ICs
  params.addRequiredParam<UserObjectName>("ics", "parameters for ics.");

  return params;
}

EntropyViscosityCoefficient::EntropyViscosityCoefficient(const InputParameters & parameters) :
    Material(parameters),
    // Boolean
    _is_dmsl_form(getParam<bool>("is_dimensional_form")),
    // Coupled variables:
    _rho(coupledValue("rho")),
    _rho_old(coupledValueOld("rho")),
    _rho_older(coupledValueOlder("rho")),
    _grad_rho(coupledGradient("rho")),
    _rhou(coupledValue("rhou")),
    _eps(isCoupled("epsilon") ?  coupledValue("epsilon") : _zero),
    _eps_old(isCoupled("epsilon") ?  coupledValueOld("epsilon") : _zero),
    _eps_older(isCoupled("epsilon") ?  coupledValueOlder("epsilon") : _zero),
    _grad_eps(isCoupled("epsilon") ? coupledGradient("epsilon") : _grad_zero),
    // Coupled aux variables
    _press(coupledValue("pressure")),
    _press_old(coupledValueOld("pressure")),
    _press_older(coupledValueOlder("pressure")),
    _grad_press(coupledGradient("pressure")),
    // Jump values:
    _jump_press(isCoupled("jump_press") ? coupledValue("jump_press") : _zero),
    _jump_dens(isCoupled("jump_dens") ? coupledValue("jump_dens") : _zero),
    // Declare material properties for viscosity coefficients.
    _kappa(declareProperty<Real>("kappa")),
    _kappa_max(declareProperty<Real>("kappa_max")),
    // Parameters
    _Cjump(getParam<double>("Cjump")),
    _is_first_order_viscosity(getParam<bool>("is_first_order_viscosity")),
    _use_jumps(getParam<bool>("use_jumps")),
    // UserObject:
    _eos(getUserObject<EquationOfState>("eos")),
    // Userobject computing the ICs
    _ics(getUserObject<InputFileSpecifiedICsRadHydro>("ics")),
    // Non-dimensional number Po
    _Po(_is_dmsl_form ? 1. : _ics.P())
{
  if (_Cjump < 0.)
    mooseError(this->name() << ": the coefficient Cjump has to be positive.");
}

void
EntropyViscosityCoefficient::computeQpProperties()
{
  // Cell size
  Real h_cell = std::pow(_current_elem->volume(),1./_mesh.dimension());

  // Compute first order viscosity:
  Real vel = _rhou[_qp]/_rho[_qp];
  Real sp = std::sqrt(_eos.c2_from_p_rho(_rho[_qp], _press[_qp], _Po*_eps[_qp]));
  if (sp < 0)
    sp = std::sqrt(_eos.c2_from_p_rho(_rho_old[_qp], _press_old[_qp], _Po*_eps_old[_qp]));
  _kappa_max[_qp] = 0.5*h_cell*(std::fabs(vel) + sp);

  // Weights for BDF2
  Real w0 = _t_step > 2 ? (2.*_dt+_dt_old)/(_dt*(_dt+_dt_old)) : 1. / _dt;
  Real w1 = _t_step > 2 ? -(_dt+_dt_old)/(_dt*_dt_old) : -1. / _dt;
  Real w2 = _t_step > 2 ? _dt/(_dt_old*(_dt+_dt_old)) : 0.;

  // Compute the jump:
  Real jump = _jump_press[_qp] + sp*sp*_jump_dens[_qp];
  Real grad = std::fabs(_grad_press[_qp](0)) + sp*sp*std::fabs(_grad_rho[_qp](0));
  Real jump_value = _use_jumps ? jump : grad;
  jump_value *= _Cjump*std::fabs(vel);

  // Compute the pressure contribution to the residual:
  Real residual=w0*_press[_qp]+w1*_press_old[_qp]+w2*_press_older[_qp];
  residual+=vel*_grad_press[_qp](0);

  // Compute the radiation contribution to the residual
  residual+=(w0*_eps[_qp]+w1*_eps_old[_qp]+w2*_eps_older[_qp])/3.;
  residual+=vel*_grad_eps[_qp](0)/3.;

  // Compute radiation contribution to the residual:
  residual-=sp*sp*(w0*_rho[_qp]+w1*_rho_old[_qp]+w2*_rho_older[_qp]);
  residual-=sp*sp*vel*_grad_rho[_qp](0);

  // Compute norm:
//  Real norm = std::min(_rho[_qp]*sp*sp, _press[_qp]);
  Real norm = 0.5*_rho[_qp]*sp*sp;

  // Compute high-order viscosity coefficient:
  Real kappa_e = _t_step == 1 ? _kappa_max[_qp] : h_cell*h_cell*(std::max(std::fabs(residual), jump_value)) / norm;

  // Get the value of the viscosity coefficients:
  _kappa[_qp] = _is_first_order_viscosity ? _kappa_max[_qp] : std::min( _kappa_max[_qp], kappa_e );
}