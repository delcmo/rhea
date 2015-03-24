#include "TimeStepCFL.h"

template<>
InputParameters validParams<TimeStepCFL>()
{
  InputParameters params = validParams<ElementPostprocessor>();

  // Coupled variables
  params.addRequiredCoupledVar("rho", "material density");
  params.addRequiredCoupledVar("rhou", "material momentum");
  params.addRequiredCoupledVar("rhoE", "material energy");
  params.addCoupledVar("radiation", "radiation");
  // Equation of state
  params.addRequiredParam<UserObjectName>("eos", "Equation of state");  
  // Parameter
  params.addParam<Real>("cfl", 0.8, "CFL number to supply by the user");

  return params;
}

TimeStepCFL::TimeStepCFL(const std::string & name, InputParameters parameters) :
    ElementPostprocessor(name, parameters),
    // Coupled variables
    _rho(coupledValue("rho")),
    _rhou(coupledValue("rhou")),
    _rhoE(coupledValue("rhoE")),
    _epsilon(isCoupled("radiation") ? coupledValue("radiation") : _zero),
    // Equation of state:
    _eos(getUserObject<EquationOfState>("eos")),
    // Parameters
    _cfl(getParam<Real>("cfl")),
    _value(0.)
{
}

TimeStepCFL::~TimeStepCFL()
{
}

void
TimeStepCFL::initialize()
{
  _value = std::numeric_limits<Real>::max();
}

void
TimeStepCFL::execute()
{
  // Compute cell size
  Real h_cell = std::pow(_current_elem->volume(), 1./_mesh.dimension());

  // Loop over quadrature points
  for (unsigned qp = 0; qp < _qrule->n_points(); ++qp)
  {
    // Compute local max eigenvalue
    Real vel = _rhou[qp] / _rho[qp];
    Real eigen = std::fabs(vel)+std::sqrt(_eos.c2(_rho[qp], vel, _rhoE[qp], _epsilon[qp]));
    Real dt = _cfl * h_cell / eigen;

    // Compute the local time step
    _value = std::min(_value, dt);
  }
}

Real
TimeStepCFL::getValue()
{
  _communicator.min(_value);
  return _value;
}

void
TimeStepCFL::threadJoin(const UserObject & uo)
{
  const TimeStepCFL & pps = dynamic_cast<const TimeStepCFL &>(uo);
  _value = std::min(_value, pps._value);
}