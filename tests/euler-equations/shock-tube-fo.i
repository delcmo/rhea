#
#####################################################
# Define some global parameters used in the blocks. #
#####################################################
#
[GlobalParams]
###### Other parameters #######
order = FIRST
isRadiation = false

###### Constans #######
speed_of_light = 1.
a = 0.
#sigma_a0 = '0. 0. 0.'
#sigma_t0 = '1. 0. 0.'
P = 1.
SIGMA_A = 0.
K = 1.
rho_hat_0 = 1.
T_hat_0 = 1.

###### Initial Conditions #######
rho_init_left = 3.
rho_init_right = 1.
vel_init_left = 0.
vel_init_right = 0.
temp_init_left = 1.
temp_init_right = 1.
eps_init_left = 0.
eps_init_right = 0.
#p_init_left = 3.
#p_init_right = 1.
membrane = 0.5
[]

#############################################################################
#                          USER OBJECTS                                     #
#############################################################################
# Define the user object class that store the EOS parameters.               #
#############################################################################

[UserObjects]
  [./eos]
    type = IdealGasEquationOfState
  	gamma = 1.4
  	Cv = 2.5
  [../]

  [./JumpGradPress]
    type = JumpGradientInterface
    variable = pressure
    jump_name = jump_grad_press
  [../]

  [./JumpGradDens]
    type = JumpGradientInterface
    variable = rho
    jump_name = jump_grad_dens
  [../]
[]

###### Mesh #######
[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 400
  xmin = 0
  xmax = 1.
  block_id = '0'
[]

#############################################################################
#                             VARIABLES                                     #
#############################################################################
# Define the variables we want to solve for: l=liquid phase,  g=vapor phase.#
#############################################################################
[Variables]
  [./rho]
    family = LAGRANGE
    scaling = 1e+0
    [./InitialCondition]
          type = InitialConditions
          eos = eos
    [../]
  [../]

  [./rhou]
    family = LAGRANGE
    scaling = 1e+0
    [./InitialCondition]
          type = InitialConditions
          eos = eos
    [../]
  [../]

  [./rhoE]
    family = LAGRANGE
    scaling = 1e+0
    [./InitialCondition]
          type = InitialConditions
          eos = eos
    [../]
  [../]
[]

############################################################################################################
#                                            KERNELS                                                       #
############################################################################################################
# Define the kernels for time dependent, convection and viscosity terms. Same index as for variable block. #
############################################################################################################

[Kernels]

  [./MassTime]
    type = TimeDerivative
    variable = rho
  [../]

  [./MomTime]
    type = TimeDerivative
    variable = rhou
  [../]

  [./EnerTime]
    type = TimeDerivative
    variable = rhoE
  [../]

  [./MassHyperbolic]
    type = RheaMass
    variable = rho
    rhou = rhou
  [../]

  [./MomHyperbloic]
    type = RheaMomentum
    variable = rhou
    rho = rho
    rhoE = rhoE
    eos = eos
  [../]

  [./EnergyHyperbolic]
    type = RheaEnergy
    variable = rhoE
    rho = rho
    rhou = rhou
    eos = eos
  [../]

  [./MassVisc]
    type = RheaArtificialVisc
    variable = rho
    equation_name = continuity
  [../]

  [./MomVisc]
    type = RheaArtificialVisc
    variable = rhou
    equation_name = x_momentum
  [../]

  [./EnergyVisc]
    type = RheaArtificialVisc
    variable = rhoE
    equation_name = energy
  [../]
[]

##############################################################################################
#                                       AUXILARY VARIABLES                                   #
##############################################################################################
# Define the auxilary variables                                                              #
##############################################################################################
[AuxVariables]
   [./pressure]
      family = LAGRANGE
   [../]

   [./temperature]
    family = LAGRANGE
   [../]

  [./mach_number]
    family = LAGRANGE
  [../]

  [./jump_grad_press]
    family = MONOMIAL
    order = CONSTANT
  [../]

  [./jump_grad_dens]
    family = MONOMIAL
    order = CONSTANT
  [../]

  [./kappa_max]
    family = MONOMIAL
    order = CONSTANT
  [../]

  [./kappa]
    family = MONOMIAL
    order = CONSTANT
  [../]
[]

##############################################################################################
#                                       AUXILARY KERNELS                                     #
##############################################################################################
# Define the auxilary kernels for liquid and gas phases. Same index as for variable block.   #
##############################################################################################
[AuxKernels]
  [./PressureAK]
    type = PressureAux
    variable = pressure
    rho = rho
    rhou = rhou
    rhoE = rhoE
    eos = eos
  [../]

  [./TemperatureAK]
    type = TemperatureAux
    variable = temperature
    rho = rho
    rhou = rhou
    rhoE = rhoE
    eos = eos
  [../]

  [./MachNumberAK]
    type = MachNumberAux
    variable = mach_number
    rho = rho
    rhou = rhou
    rhoE = rhoE
    eos = eos
  [../]

  [./KappaMaxAK]
    type = MaterialRealAux
    variable = kappa_max
    property = kappa_max
  [../]

  [./KappaAK]
    type = MaterialRealAux
    variable = kappa
    property = kappa
  [../]
[]

##############################################################################################
#                                       MATERIALS                                            #
##############################################################################################
# Define functions that are used in the kernels and aux. kernels.                            #
##############################################################################################
[Materials]
  [./EntViscCoeff]
    type = EntropyViscosityCoefficient
    block = '0'
    rho = rho
    rhou = rhou
    pressure = pressure
    jump_press = jump_grad_press
    jump_dens = jump_grad_dens
    Cjump = 1.
    is_first_order_viscosity = true
    eos = eos
  [../]

  [./PhysicalPropertyMaterial]
    type = PhysicalPropertyMaterial
    block = '0'
    rho = rho
    pressure = pressure
    eos = eos
  [../]
[]
##############################################################################################
#                               BOUNDARY CONDITIONS                                          #
##############################################################################################
# Define the functions computing the inflow and outflow boundary conditions.                 #
##############################################################################################
[BCs]
  [./MassLeft]
    type = DirichletBC
    variable = rho
    value = 3.
    boundary = 'left'
  [../]

  [./MassRight]
    type = DirichletBC
    variable = rho
    value = 1.
    boundary = 'right'
  [../]

  [./MomentumLeft]
    type = DirichletBC
    variable = rhou
    value = 0.
    boundary = 'left'
  [../]

  [./MomentumRight]
    type = DirichletBC
    variable = rhou
    value = 0.
    boundary = 'right'
  [../]

  [./EnergyLeft]
    type = DirichletBC
    variable = rhoE
    value = 7.5
    boundary = 'left'
  [../]

  [./EnergyRight]
    type = DirichletBC
    variable = rhoE
    value = 2.5
    boundary = 'right'
  [../]
[]
##############################################################################################
#                                  PRECONDITIONER                                            #
##############################################################################################
# Define the functions computing the inflow and outflow boundary conditions.                 #
##############################################################################################

[Preconditioning]
  active = 'FDP_Newton'
  [./FDP_Newton]
    type = FDP
    full = true
    solve_type = 'PJFNK'
    #petsc_options = '-snes_mf_operator -snes_ksp_ew'
    #petsc_options_iname = '-mat_fd_coloring_err  -mat_fd_type  -mat_mffd_type'
    #petsc_options_value = '1.e-12       ds             ds'
  [../]
[]

##############################################################################################
#                                     EXECUTIONER                                            #
##############################################################################################
# Define the functions computing the inflow and outflow boundary conditions.                 #
##############################################################################################

[Executioner]
  type = Transient
  scheme = 'bdf2'
  end_time = 1.
  dt = 1.e-4
  dtmin = 1e-9
  l_tol = 1e-8
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-7
  l_max_its = 50
  nl_max_its = 50
  num_steps = 50
  [./TimeStepper]
    type = FunctionDT
    time_t =  '0.     3.e-2   1.e-1  1.'
    time_dt = '1.e-3  1.e-3   1.e-3  1.e-3'
  [../]
[]

##############################################################################################
#                                        OUTPUT                                              #
##############################################################################################
# Define the functions computing the inflow and outflow boundary conditions.                 #
##############################################################################################

[Outputs]
  output_initial = true
  interval = 1
  exodus = true
[]

##############################################################################################
#                                        DEBUG                                               #
##############################################################################################
# Debug                 #
##############################################################################################

#[Debug]
#  show_var_residual_norms = true
#[]
