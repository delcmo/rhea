#
#####################################################
# Define some global parameters used in the blocks. #
#####################################################
#
[GlobalParams]
###### Other parameters #######
order = FIRST

###### Constans #######
speed_of_light = 1.
a = 1.
cross_section_name = temp_dpt_cs
is_diffusion = false
sigma_a0 = '0. 0. 0.'
sigma_t0 = '2.e+003 0. 0.'
rho_hat_0 = 1.
T_hat_0 = 1.
K = 1.
SIGMA_A = 1.
P = 1.
Mach_inlet = 1.

###### Initial Conditions #######
rho_init_left = 1.
rho_init_right = 1.
vel_init_left = 2.
vel_init_right = 2.
temp_init_left = 1.
temp_init_right = 1.
eps_init_left = 3.
eps_init_right = 1.
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

  [./ics]
    type = ComputeICsRadHydro
    eos = eos
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
  nx = 200
  xmin = 0.
  xmax = 1.
  block_id = '0'
[]

#############################################################################
#                             VARIABLES                                     #
#############################################################################
# Define the variables we want to solve for: l=liquid phase,  g=vapor phase.#
#############################################################################
[Variables]
  [./epsilon]
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
  [./EpsilonTime]
    type = TimeDerivative
    variable = epsilon
  [../]

  [./RadiationHyperbolic]
    type = RheaRadiation
    variable = epsilon
    rho = rho
    rhou = rhou
    rhoE = rhoE
    eos = eos
    ics = ics
  [../]

  [./RadiationVisc]
    type = RheaArtificialVisc
    variable = epsilon
    equation_name = radiation
  [../]
[]

##############################################################################################
#                                       AUXILARY VARIABLES                                   #
##############################################################################################
# Define the auxilary variables                                                              #
##############################################################################################
[AuxVariables]
   [./rho]
      family = LAGRANGE
      [./InitialCondition]
        type = ConstantIC
        value = 1.
      [../]
   [../]

   [./rhou]
    family = LAGRANGE
      [./InitialCondition]
        type = ConstantIC
        value = 2.
      [../]
   [../]
   
   [./rhoE]
    family = LAGRANGE
      [./InitialCondition]
        type = ConstantIC
        value = 1.
      [../]
   [../]
   
   [./pressure]
    family = LAGRANGE
      [./InitialCondition]
        type = ConstantIC
        value = 1.
      [../]
   [../]

  [./diffusion]
    family = MONOMIAL
    order = CONSTANT
  [../]

  [./jump_grad_dens]
    family = MONOMIAL
    order = CONSTANT
  [../]

  [./jump_grad_press]
    family = MONOMIAL
    order = CONSTANT
  [../]

  [./kappa]
    family = MONOMIAL
    order = CONSTANT
  [../]
  
  [./kappa_max]
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
  [./DiffusionAK]
    type = MaterialRealAux
    variable = diffusion
    property = diffusion
  [../]

  [./KappaAK]
    type = MaterialRealAux
    variable = kappa
    property = kappa
  [../]
  
  [./KappaMaxAK]
    type = MaterialRealAux
    variable = kappa_max
    property = kappa_max
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
    epsilon = epsilon
    pressure = pressure
    jump_press = jump_grad_press
    jump_dens = jump_grad_dens
    Cjump = 15.
    is_first_order_viscosity = false
    eos = eos
  [../]

  [./PhysicalPropertyMaterial]
    type = PhysicalPropertyMaterial
    block = '0'
    rho = rho
    pressure = pressure
    eos = eos
    ics = ics    
  [../]
[]
##############################################################################################
#                               BOUNDARY CONDITIONS                                          #
##############################################################################################
# Define the functions computing the inflow and outflow boundary conditions.                 #
##############################################################################################
[BCs]
  [./RadiationRight]
    type = DirichletBC
    variable = epsilon
    value = 3.
    boundary = 'left'
  [../]
  
  [./RadiationLeft]
    type = DirichletBC
    variable = epsilon
    value = 1.
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
  l_max_its = 10
  nl_max_its = 10
  num_steps = 25
  [./TimeStepper]
    type = FunctionDT
    time_t =  '0.     1.'
    time_dt = '1.e-3  1.e-3'
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
