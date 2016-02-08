#include "RheaApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

// Kernels
#include "RheaMass.h"
#include "RheaMomentum.h"
#include "RheaEnergy.h"
#include "RheaRadiation.h"
#include "RheaNonDimensionalMomentum.h"
#include "RheaNonDimensionalEnergy.h"
#include "RheaNonDimensionalRadiation.h"
#include "RheaArtificialVisc.h"
#include "RheaForcingFunction.h"
#include "RheaForcingFunctionStream.h"

// Auxkernels
#include "PressureAux.h"
#include "TemperatureAux.h"
#include "InternalEnergyAux.h"
#include "VelocityAux.h"
#include "MachNumberAux.h"
#include "RadTempAux.h"

// BCs
#include "RheaBCs.h"

// ICs
#include "InitialConditions.h"
#include "RheaIC.h"

// Userobjects
#include "EquationOfState.h"
#include "IdealGasEquationOfState.h"
#include "ComputeICsRadHydro.h"
#include "InputFileSpecifiedICsRadHydro.h"
#include "JumpGradientInterface.h"

//Materials
#include "PhysicalPropertyMaterial.h"
#include "EntropyViscosityCoefficient.h"
#include "NonDimensionalPhysicalPropertyMaterial.h"
#include "NonDimensionalEntropyViscosityCoefficient.h"

// Postprocessors
#include "TimeStepCFL.h"

template<>
InputParameters validParams<RheaApp>()
{
  InputParameters params = validParams<MooseApp>();

  params.set<bool>("use_legacy_uo_initialization") = false;
  params.set<bool>("use_legacy_uo_aux_computation") = false;
  params.set<bool>("use_legacy_output_syntax") = false;
  
  return params;
}

RheaApp::RheaApp(InputParameters parameters) :
    MooseApp(parameters)
{
  Moose::registerObjects(_factory);
  ModulesApp::registerObjects(_factory);
  RheaApp::registerObjects(_factory);

  Moose::associateSyntax(_syntax, _action_factory);
  ModulesApp::associateSyntax(_syntax, _action_factory);
  RheaApp::associateSyntax(_syntax, _action_factory);
}

RheaApp::~RheaApp()
{
}

extern "C" void RheaApp__registerApps() { RheaApp::registerApps(); }
void
RheaApp::registerApps()
{
  registerApp(RheaApp);
}

extern "C" void RheaApp__registerObjects(Factory & factory) { RheaApp::registerObjects(factory); }
void
RheaApp::registerObjects(Factory & factory)
{
  // Kernels
  registerKernel(RheaMass);
  registerKernel(RheaMomentum);
  registerKernel(RheaEnergy);
  registerKernel(RheaRadiation);
  registerKernel(RheaNonDimensionalMomentum);
  registerKernel(RheaNonDimensionalEnergy);
  registerKernel(RheaNonDimensionalRadiation);
  registerKernel(RheaArtificialVisc);
  registerKernel(RheaForcingFunction);
  registerKernel(RheaForcingFunctionStream);

  // Auxkernels
  registerAux(PressureAux);
  registerAux(TemperatureAux);
  registerAux(InternalEnergyAux);
  registerAux(VelocityAux);
  registerAux(MachNumberAux);
  registerAux(RadTempAux);

  // BCs
  registerBoundaryCondition(RheaBCs);

  // ICs
  registerInitialCondition(InitialConditions);
  registerInitialCondition(RheaIC);

  // Userobjects
  registerUserObject(EquationOfState);
  registerUserObject(IdealGasEquationOfState);
  registerUserObject(ComputeICsRadHydro);
  registerUserObject(InputFileSpecifiedICsRadHydro);
  registerUserObject(JumpGradientInterface);

  // Materials
  registerMaterial(PhysicalPropertyMaterial);
  registerMaterial(EntropyViscosityCoefficient);
  registerMaterial(NonDimensionalPhysicalPropertyMaterial);
  registerMaterial(NonDimensionalEntropyViscosityCoefficient);

  // Postprocessors
  registerPostprocessor(TimeStepCFL);
}

extern "C" void RheaApp__associateSyntax(Syntax & syntax, ActionFactory & action_factory) { RheaApp::associateSyntax(syntax, action_factory); }
void
RheaApp::associateSyntax(Syntax & /*syntax*/, ActionFactory & /*action_factory*/)
{
}