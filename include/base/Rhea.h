#ifndef RHEA_H
#define RHEA_H

class Factory;
class ActionFactory;
class Syntax;

namespace Rhea
{
  /**
   * Register this application and any it depends on.
   */
  void registerApps();

  /**
   * Registers all Kernels and BCs, etc.
   */
  void registerObjects(Factory & factory);

  /**
   * Register syntax and Action objects.
   */
  void associateSyntax(Syntax & syntax, ActionFactory & action_factory);
}

#endif //RHEA_H
