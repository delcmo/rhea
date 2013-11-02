#
#####################################################
# Define some global parameters used in the blocks. #
#####################################################
#
[GlobalParams]
###### Other parameters #######
order = FIRST
viscosity_name = ENTROPY
isRadiation = false
Ce = 1.

###### Constans #######
speed_of_light = 299.792
a = 1.372e-2
sigma_a0 = 0.
sigma_t0 = 8.5314410158161809e+002

###### Initial Conditions #######
rho_init_left = 15.
rho_init_right = 1.
vel_init_left = 0
vel_init_right = 0
temp_init_left = 1
temp_init_right = 1
eps_init_left = 0.
#p_init_left = 1. # 3.
#p_init_right = 1
eps_init_right = 0.
membrane = 0.7
[]

#############################################################################
#                          USER OBJECTS                                     #
#############################################################################
# Define the user object class that store the EOS parameters.               #
#############################################################################

[UserObjects]
  [./eos]
    type = EquationOfState
  	gamma = 1.4
  	Pinf = 0.
  	q = 0.
  	Cv = 2.5
  	q_prime = 0.
  [../]

  [./JumpGradPress]
    type = JumpGradientInterface
    variable = pressure
    jump_name = jump_grad_press
    execute_on = timestep_begin
  [../]
[]

###### Mesh #######
[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 100
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
    type = RheaTimeDerivative
    variable = rho
  [../]

  [./MomTime]
    type = RheaTimeDerivative
    variable = rhou
  [../]

  [./EnerTime]
    type = RheaTimeDerivative
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
    velocity = velocity
    pressure = pressure
  [../]

  [./EnergyHyperbolic]
    type = RheaEnergy
    variable = rhoE
    velocity = velocity
    temperature = temperature
    pressure = pressure
  [../]

  [./MassVisc]
    type = RheaArtificialVisc
    variable = rho
    equation_name = CONTINUITY
    density = rho
    internal_energy = internal_energy
    velocity = velocity
  [../]

  [./MomVisc]
    type = RheaArtificialVisc
    variable = rhou
    equation_name = MOMENTUM
    density = rho
    internal_energy = internal_energy
    velocity = velocity
  [../]

  [./EnergyVisc]
    type = RheaArtificialVisc
    variable = rhoE
    equation_name = ENERGY
    density = rho
    internal_energy = internal_energy
    velocity = velocity
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
   [../]

   [./internal_energy]
      family = LAGRANGE
   [../]

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

  [./mu_max]
    family = MONOMIAL
    order = CONSTANT
  [../]

  [./kappa_max]
    family = MONOMIAL
    order = CONSTANT
  [../]

  [./mu]
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
 [./VelocityAK]
    type = VelocityAux
    variable = velocity
    rho = rho
    rhou = rhou
  [../]

  [./IntEnerAK]
    type = InternalEnergyAux
    variable = internal_energy
    rho = rho
    rhou = rhou
    rhoE = rhoE
  [../]

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
    pressure = pressure
    rho = rho
    eos = eos
  [../]

  [./RadTempAK]
    type = RadTempAux
    variable = rad_temp
    temperature = temperature
  [../]

  [./MachNumberAK]
    type = MachNumberAux
    variable = mach_number
    rho = rho
    rhou = rhou
    pressure = pressure
    eos = eos
  [../]

  [./MuMaxAK]
    type = MaterialRealAux
    variable = mu_max
    property = mu_max
  [../]

  [./KappaMaxAK]
    type = MaterialRealAux
    variable = kappa_max
    property = kappa_max
  [../]

  [./MuAK]
    type = MaterialRealAux
    variable = mu
    property = mu
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
  [./ViscCoeff]
    type = ComputeMaterials
    block = '0'
    velocity = velocity
    pressure = pressure
    density = rho
    jump = jump_grad_press
    velocity_PPS_name = MaxVelocity
    eos = eos
  [../]
[]

##############################################################################################
#                                     PPS                                                    #
##############################################################################################
# Define functions that are used in the kernels and aux. kernels.                            #
##############################################################################################
[Postprocessors]

[./MaxVelocity]
    type = NodalMaxValue
    variable = velocity
    execute_on = timestep_begin
[../]

[]
##############################################################################################
#                               BOUNDARY CONDITIONS                                          #
##############################################################################################
# Define the functions computing the inflow and outflow boundary conditions.                 #
##############################################################################################
[BCs]
  [./MassLeft]
   type = RheaBCs
    variable = rho
    equation_name = CONTINUITY
    velocity = velocity
    density = rho
    pressure = pressure
    eos = eos
    p_bc = 15
    T_bc = 1
    v_bc = 0
#    type = DirichletBC
#    value = 3.
    boundary = 'left'
  [../]

  [./MassRight]
   type = RheaBCs
    variable = rho
    equation_name = CONTINUITY
   velocity = velocity
    density = rho
    pressure = pressure
    eos = eos
    p_bc = 1
    T_bc = 1
#    type = DirichletBC
#    value = 1.
    boundary = 'right'
  [../]

  [./MomentumLeft]
    type = RheaBCs
    variable = rhou
    equation_name = MOMENTUM
    velocity = velocity
    density = rho
    pressure = pressure
    eos = eos
    p_bc = 15
    T_bc = 1
    v_bc = 0
#    type = DirichletBC
#    value = 0.
    boundary = 'left'
  [../]

  [./MomentumRight]
    type = RheaBCs
    variable = rhou
    equation_name = MOMENTUM
    velocity = velocity
    density = rho
    pressure = pressure
    eos = eos
    p_bc = 1
    T_bc = 1
#    type = DirichletBC
#    value = 0.
    boundary = 'right'
  [../]

  [./EnergyLeft]
    type = RheaBCs
    variable = rhoE
    equation_name = ENERGY
    velocity = velocity
    density = rho
    pressure = pressure
    eos = eos
    p_bc = 15
    T_bc = 1
    v_bc = 0
#    type = DirichletBC
#    value = 7.5
    boundary = 'left'
  [../]

  [./EnergyRight]
    type = RheaBCs
    variable = rhoE
    equation_name = ENERGY
    velocity = velocity
    density = rho
    pressure = pressure
    eos = eos
    p_bc = 1
    T_bc = 1
#    type = DirichletBC
#    value = 2.5
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
    petsc_options_iname = '-mat_fd_coloring_err  -mat_fd_type  -mat_mffd_type'
    petsc_options_value = '1.e-12       ds             ds'
  [../]
[]

##############################################################################################
#                                     EXECUTIONER                                            #
##############################################################################################
# Define the functions computing the inflow and outflow boundary conditions.                 #
##############################################################################################

[Executioner]
  type = Transient
  string scheme = 'bdf2'
  num_steps = 300
  #end_time = 200
  dt = 1.e-3
  dtmin = 1e-9
  l_tol = 1e-8
  nl_rel_tol = 1e-8
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
