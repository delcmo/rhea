#ifndef RHEAAPP_H
#define RHEAAPP_H

#include "MooseApp.h"

class RheaApp;

template<>
InputParameters validParams<RheaApp>();

class RheaApp : public MooseApp
{
public:
  RheaApp(InputParameters parameters);
  virtual ~RheaApp();

  static void registerApps();
  static void registerObjects(Factory & factory);
  static void associateSyntax(Syntax & syntax, ActionFactory & action_factory);
};

#endif /* RHEAAPP_H */
