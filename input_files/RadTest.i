#
#####################################################
# Define some global parameters used in the blocks. #
#####################################################
#
[GlobalParams]
###### Other parameters #######
order = FIRST

###### Constans #######
speed_of_light = 299.792
a = 1.372e-2
sigma_a0 = 3.9071164263502112e+002
sigma_t0 = 8.5314410158161809e+002

###### Initial Conditions #######
rho_init_left = 1.
rho_init_right = 1.
vel_init_left = 1.2298902390050911e-001
vel_init_right = 1.2298902390050911e-001
temp_init_left = 1.0000000000000001e-001
temp_init_right = 1.0000000000000001e-001
eps_init_left = 2
p_init_left = 8.232e-03
p_init_right = 8.232e-03
eps_init_right = 2
membrane = 0.5
eps_left = 0.
eps_right = 0.
[]

#############################################################################
#                          USER OBJECTS                                     #
#############################################################################
# Define the user object class that store the EOS parameters.               #
#############################################################################

[UserObjects]
  [./eos]
    type = EquationOfState
  	gamma = 1.6666667
  	Pinf = 0.
  	q = 0.
  	Cv = 1.2348000000000001e-001
  	q_prime = 0.
  [../]

[]

###### Mesh #######
[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 200
  xmin = 0
  xmax = 1
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
    scaling = 1e+2
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
    type = RheaTimeDerivative
    variable = epsilon
  [../]

  [./RadiationHyperbolic]
    type = RheaRadiation
    variable = epsilon
    velocity = velocity
    temperature = temperature
  [../]
[]

##############################################################################################
#                                       AUXILARY VARIABLES                                   #
##############################################################################################
# Define the auxilary variables                                                              #
##############################################################################################
[AuxVariables]
   [./velocity]
      family = LAGRANGE
      [./InitialCondition]
        type = ConstantIC
        value = 0.
      [../]
   [../]

  [./temperature]
    family = LAGRANGE
    [./InitialCondition]
    type = ConstantIC
    value = 0.
    [../]
  [../]

  [./pressure]
    family = LAGRANGE
    [./InitialCondition]
    type = ConstantIC
    value = 1.
    [../]
  [../]

  [./rho]
    family = LAGRANGE
    [./InitialCondition]
    type = ConstantIC
    value = 1.
    [../]
  [../]
[]

##############################################################################################
#                                       AUXILARY KERNELS                                     #
##############################################################################################
# Define the auxilary kernels for liquid and gas phases. Same index as for variable block.   #
##############################################################################################
[AuxKernels]
[]

##############################################################################################
#                                       MATERIALS                                            #
##############################################################################################
# Define functions that are used in the kernels and aux. kernels.                            #
##############################################################################################

[Materials]
  [./ViscCoeff]
    type = ComputeMaterials
    block = '0'
    velocity = velocity
    pressure = pressure
    density = rho
    epsilon = epsilon
    pressure_PPS_name = AveragePressure
    velocity_PPS_name = AverageVelocity
    eos = eos
    Ce = 1.
  [../]
[]

##############################################################################################
#                                     PPS                                                    #
##############################################################################################
# Define functions that are used in the kernels and aux. kernels.                            #
##############################################################################################
[Postprocessors]

[./AverageVelocity]
    type = ElementAverageValue
    variable = velocity
[../]

[./AveragePressure]
    type = ElementAverageValue
    variable = pressure
[../]

[]
##############################################################################################
#                               BOUNDARY CONDITIONS                                          #
##############################################################################################
# Define the functions computing the inflow and outflow boundary conditions.                 #
##############################################################################################
[BCs]
  [./RadiationLeft]
    type = RheaBCs
    variable = epsilon
    equation_name = RADIATION
    velocity = velocity
    temperature = temperature
    density = rho
    pressure = pressure
    eos = eos
    p_bc = 8.232e-03
    T_bc = 1.0000000000000001e-001
    v_bc = 1.2298902390050911e-001
    boundary = 'left'
  [../]

  [./RadiationRight]
    type = RheaBCs
    variable = epsilon
    equation_name = RADIATION
    velocity = velocity
    temperature = temperature
    density = rho
    pressure = pressure
    eos = eos
    p_bc = 8.232e-03
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
    petsc_options = '-snes_mf_operator -snes_ksp_ew'
    petsc_options_iname = '-mat_fd_coloring_err'
    petsc_options_value = '1.e-13'#'1.e-12'
    #petsc_options = '-snes_mf_operator -ksp_converged_reason -ksp_monitor -snes_ksp_ew' 
    #petsc_options_iname = '-pc_type'
    #petsc_options_value = 'lu'
  [../]
[]

##############################################################################################
#                                     EXECUTIONER                                            #
##############################################################################################
# Define the functions computing the inflow and outflow boundary conditions.                 #
##############################################################################################

[Executioner]
  type = Transient   # Here we use the Transient Executioner
  #scheme = 'explicit-euler'
  string scheme = 'bdf2'
  #petsc_options = '-snes'
  #petsc_options_iname = '-pc_type'
  #petsc_options_value = 'lu'
  num_steps = 10
  #end_time = 200
  dt = 1e-1
  dtmin = 1e-9
  l_tol = 1e-8
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-6
  l_max_its = 50
  nl_max_its = 50
#  [./TimeStepper]
#    type = FunctionDT
#    time_t =  '0.     2.6e-2  5.e-1  0.56'
#    time_dt = '1e-4   1e-4    1e-3    1e-3'
#  [../]
[]

##############################################################################################
#                                        OUTPUT                                              #
##############################################################################################
# Define the functions computing the inflow and outflow boundary conditions.                 #
##############################################################################################

[Output]
  output_initial = true
  interval = 1
  exodus = true
  postprocessor_screen = false
  perf_log = true
[]

##############################################################################################
#                                        DEBUG                                               #
##############################################################################################
# Debug                 #
##############################################################################################

#[Debug]
#  show_var_residual_norms = true
#[]
