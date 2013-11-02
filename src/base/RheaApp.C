#include "RheaApp.h"
#include "Rhea.h"
#include "Moose.h"
//#include "Elk.h"

template<>
InputParameters validParams<RheaApp>()
{
  InputParameters params = validParams<MooseApp>();
  return params;
}

RheaApp::RheaApp(const std::string & name, InputParameters parameters) :
    MooseApp(name, parameters)
{
  srand(libMesh::processor_id());
  
  Moose::registerObjects(_factory);
  //Elk::registerObjects(_factory);
  Rhea::registerObjects(_factory);

  Moose::associateSyntax(_syntax, _action_factory);
  //Elk::associateSyntax(_syntax, _action_factory);
  Rhea::associateSyntax(_syntax, _action_factory);
}

RheaApp::~RheaApp()
{
}

