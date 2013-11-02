#include "Rhea.h"
#include "RheaApp.h"
#include "Moose.h"
#include "Factory.h"
#include "AppFactory.h"

// Kernels
#include "RheaTimeDerivative.h"
#include "RheaMass.h"
#include "RheaMomentum.h"
#include "RheaEnergy.h"
#include "RheaRadiation.h"
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

// Userobjects
#include "EquationOfState.h"
#include "JumpGradientInterface.h"

//Materials
#include "ComputeMaterials.h"

namespace Rhea
{
  void registerApps()
  {
    registerApp(RheaApp);
  }

  void registerObjects(Factory & factory)
  {
      // Kernels
      registerKernel(RheaTimeDerivative);
      registerKernel(RheaMass);
      registerKernel(RheaMomentum);
      registerKernel(RheaEnergy);
      registerKernel(RheaRadiation);
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
      
      // Userobjects
      registerUserObject(EquationOfState);
      registerUserObject(JumpGradientInterface);
      
      // Materials
      registerMaterial(ComputeMaterials);
  }

  void associateSyntax(Syntax & syntax, ActionFactory & action_factory)
  {
  }
}
