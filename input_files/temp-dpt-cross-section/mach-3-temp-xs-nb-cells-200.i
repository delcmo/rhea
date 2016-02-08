#
#####################################################
# Define some global parameters used in the blocks. #
#####################################################
#
[GlobalParams]
###### Other parameters #######
order = FIRST
Cjump = 2.
is_first_order_viscosity = false
use_jumps = false
cfl = 10

###### Constants #######
cross_section_name = temp_dpt_opacity
speed_of_light = 2.99792e+2
a = 1.372e-2
sigma_a0 = '7216.875 0. 3.5'
sigma_t0 = '7216.875 0. 3.5'

###### Initial Conditions #######
Mach_inlet = 3.
rho_hat_0 = 1.
T_hat_0 = 0.12156013625
P = 1.e-4
K = 1.
SIGMA_A = 1.e6
membrane = 0.
[]

#############################################################################
#                          USER OBJECTS                                     #
#############################################################################
# Define the user object class that store the EOS parameters.               #
#############################################################################

[UserObjects]
  [./eos]
    type = IdealGasEquationOfState
  	gamma = 1.6666667
  	Cv = 0.221804 # 1.2348000000000001e-001
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
  nx =200
  xmin = -5.e-02 # -2.e-2
  xmax = 5.e-02 # 2.e-2
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
    scaling = 1e+2
    [./InitialCondition]
      type = RheaIC
      eos = eos
      ics = ics
    [../]
  [../]

  [./rhou]
    family = LAGRANGE
    scaling = 1e+2
    [./InitialCondition]
      type = RheaIC
      eos = eos
      ics = ics
    [../]
  [../]

  [./rhoE]
    family = LAGRANGE
    scaling = 1e+2
    [./InitialCondition]
      type = RheaIC
      eos = eos
      ics = ics
    [../]
  [../]

  [./epsilon]
    family = LAGRANGE
    scaling = 1e+2
    [./InitialCondition]
      type = RheaIC
      eos = eos
      ics = ics      
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

  [./EpsilonTime]
    type = TimeDerivative
    variable = epsilon
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
    radiation = epsilon
    eos = eos
    ics = ics
  [../]

  [./EnergyHyperbolic]
    type = RheaEnergy
    variable = rhoE
    rho = rho
    rhou = rhou
    radiation = epsilon
    eos = eos
    ics = ics
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
   [./pressure]
      family = LAGRANGE
   [../]

   [./temperature]
    family = LAGRANGE
   [../]

  [./rad_temp]
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

  [./diffusion]
    family = MONOMIAL
    order = CONSTANT
  [../]

  [./sigma_a]
    family = MONOMIAL
    order = CONSTANT
  [../]

  [./sigma_t]
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

  [./RadTempAK]
    type = RadTempAux
    variable = rad_temp
    radiation = epsilon
  [../]

  [./MachNumberAK]
    type = MachNumberAux
    variable = mach_number
    rho = rho
    rhou = rhou
    rhoE = rhoE
    epsilon = epsilon
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

  [./DiffusionAK]
    type = MaterialRealAux
    variable = diffusion
    property = diffusion
  [../]

  [./SigmaA_AK]
  type = MaterialRealAux
    variable = sigma_a
    property = sigma_a
  [../]

  [./SigmaT_AK]
  type = MaterialRealAux
    variable = sigma_t
    property = sigma_t
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
    eos = eos
    ics = ics
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
  [./MassRight]
    type = RheaBCs
    variable = rho
    equation_name = continuity
    rho = rho
    rhou = rhou
    rhoE = rhoE
    epsilon = epsilon
    pressure = pressure
    eos = eos
    ics = ics
    boundary = 'right left'
  [../]

  [./MomentumRight]
    type = RheaBCs
    variable = rhou
    equation_name = x_momentum
    rho = rho
    rhou = rhou
    rhoE = rhoE
    epsilon = epsilon
    pressure = pressure
    eos = eos
    ics = ics
    boundary = 'right left'
  [../]

  [./EnergyRight]
    type = RheaBCs
    variable = rhoE
    equation_name = energy
    rho = rho
    rhou = rhou
    rhoE = rhoE
    epsilon = epsilon
    pressure = pressure
    eos = eos
    ics = ics
    boundary = 'right left'
  [../]
  
  [./RadiationLeft]
    type = DirichletBC
    variable = epsilon
    value = 2.995841442e-06
    boundary = 'left'
  [../]

  [./RadiationRight]
    type = DirichletBC
    variable = epsilon
    value = 0.0005387047241
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
#                                  POSTPROCESSORS                                            #
##############################################################################################
# Define the functions computing the inflow and outflow boundary conditions.                 #
##############################################################################################

[Postprocessors]
  [./dt]
    type = TimeStepCFL
    rho = rho
    rhou = rhou
    rhoE = rhoE
    radiation = epsilon
    eos = eos
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
  end_time = 10.
  dt = 1.e-4
  dtmin = 1e-9
  l_tol = 1e-8
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-7
  l_max_its = 50
  nl_max_its = 50
#  trans_ss_check = true
#  ss_check_tol = 1.e-12
  [./TimeStepper]
    type = PostprocessorDT
    postprocessor = dt
    dt = 1.e-4
  [../]
[]

##############################################################################################
#                                        OUTPUT                                              #
##############################################################################################
# Define the functions computing the inflow and outflow boundary conditions.                 #
##############################################################################################

[Outputs]
  file_base = mach-3-temp-xs-nb-cells-200_out
  [./console]
  type = Console
  perf_log = true
  interval = 1
  [../]
  
  [./out]
    type = Exodus
    interval = 100
    output_initial = true    
  [../]
[]

##############################################################################################
#                                        DEBUG                                               #
##############################################################################################
# Debug                 #
##############################################################################################

#[Debug]
#  show_var_residual_norms = true
#[]
