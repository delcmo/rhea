#ifndef RHEAAPP_H
#define RHEAAPP_H

#include "MooseApp.h"

class RheaApp;

template<>
InputParameters validParams<RheaApp>();

class RheaApp : public MooseApp
{
public:
  RheaApp(const std::string & name, InputParameters parameters);
  virtual ~RheaApp();
};

#endif /* RHEAAPP_H */
