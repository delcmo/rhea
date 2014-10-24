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

###### Constans #######
speed_of_light = 299.792
a = 1.372e-2
sigma_a0 = 3.9071164263502112e+002 # 002
sigma_t0 = 8.5314410158161809e+008 # 002

###### Initial Conditions #######
rho_init_left = 1.
rho_init_right = 1.0749588725963066e+000 
vel_init_left = 1.2298902390050911e-001
vel_init_right = 1.1441277153558302e-001
temp_init_left = 1.0000000000000001e-001
temp_init_right = 1.0494545175154467e-001
eps_init_left = 1.3720000000000002e-006
eps_init_right = 1.6642117992569650e-006
#p_init_left = 8.232e-03
#p_init_right = 9.28508e-03
membrane = 0.
length = 1e-2
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
  xmin = -5.27105579632e-002
  xmax = 4.62235483783162e-002
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
    scaling = 1e+4
	[./InitialCondition]
        type = InitialConditions
        eos = eos
	[../]
  [../]

  [./rhou]
    family = LAGRANGE
    scaling = 1e+4
	[./InitialCondition]
        type = InitialConditions
        eos = eos
	[../]
  [../]

  [./rhoE]
    family = LAGRANGE
    scaling = 1e+4
	[./InitialCondition]
        type = InitialConditions
        eos = eos
	[../]
  [../]

  [./epsilon]
    family = LAGRANGE
    scaling = 1e+1
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

  [./EpsilonTime]
    type = RheaTimeDerivative
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
    velocity = velocity
    pressure = pressure
    radiation = epsilon
  [../]

  [./EnergyHyperbolic]
    type = RheaEnergy
    variable = rhoE
    velocity = velocity
    temperature = temperature
    pressure = pressure
    radiation = epsilon
  [../]

  [./RadiationHyperbolic]
    type = RheaRadiation
    variable = epsilon
    velocity = velocity
    temperature = temperature
  [../]

  [./MassVisc]
    type = RheaArtificialVisc
    variable = rho
    equation_name = CONTINUITY
    density = rho
    internal_energy = internal_energy
    velocity = velocity
    radiation = epsilon
  [../]

  [./MomVisc]
    type = RheaArtificialVisc
    variable = rhou
    equation_name = MOMENTUM
    density = rho
    internal_energy = internal_energy
    velocity = velocity
    radiation = epsilon
  [../]

  [./EnergyVisc]
    type = RheaArtificialVisc
    variable = rhoE
    equation_name = ENERGY
    density = rho
    internal_energy = internal_energy
    velocity = velocity
    radiation = epsilon
  [../]

  [./RadiationVisc]
    type = RheaArtificialVisc
    variable = epsilon
    equation_name = RADIATION
    density = rho
    internal_energy = internal_energy
    velocity = velocity
    radiation = epsilon
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

  [./radiation_temperature]
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

  [./diffusion]
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
    variable = radiation_temperature
    temperature = epsilon
  [../]

  [./MachNumberAK]
    type = MachNumberAux
    variable = mach_number
    rho = rho
    rhou = rhou
    pressure = pressure
    eos = eos
  [../]

  [./DiffusionAK]
    type = MaterialRealAux
    variable = diffusion
    property = diffusion
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
  [./ViscCoeff]
    type = ComputeMaterials
    block = '0'
    velocity = velocity
    pressure = pressure
    density = rho
    epsilon = epsilon
    jump_press = jump_grad_press
    jump_dens = jump_grad_dens
    eos = eos
    Ce = 1.2
  [../]
[]

##############################################################################################
#                                     PPS                                                    #
##############################################################################################
# Define functions that are used in the kernels and aux. kernels.                            #
##############################################################################################
[Postprocessors]

#[./AverageVelocity]
#    type = ElementAverageValue
#    variable = velocity
#[../]

#[./AverageEpsilon]
#    type = ElementAverageValue
#    variable = epsilon
#[../]

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
    value = 1.
    boundary = 'left'
  [../]

  [./MassRight]
    type = RheaBCs
    variable = rho
    equation_name = CONTINUITY
    velocity = velocity
    density = rho
    pressure = pressure
    epsilon = epsilon
    eos = eos
    p_bc = 9.28508e-03
    boundary = 'right'
  [../]

  [./MomentumLeft]
    type = DirichletBC
    variable = rhou
    value = 0.122989024
    boundary = 'left'
  [../]

  [./MomentumRight]
    type = RheaBCs
    variable = rhou
    equation_name = MOMENTUM
    velocity = velocity
    density = rho
    pressure = pressure
    epsilon = epsilon
    eos = eos
    p_bc = 9.28508e-03
    boundary = 'right'
  [../]

  [./EnergyLeft]
    type = DirichletBC
    variable = rhoE
    value = 0.01991115
    boundary = 'left'
  [../]

  [./EnergyRight]
    type = RheaBCs
    variable = rhoE
    equation_name = ENERGY
    velocity = velocity
    density = rho
    pressure = pressure
    epsilon = epsilon
    eos = eos
    p_bc = 9.28508e-03
    boundary = 'right'
  [../]

  [./RadiationLeft]
#    type = DirichletBC
    type = RheaBCs
    variable = epsilon
    equation_name = RADIATION
    velocity = velocity
    density = rho
    pressure = pressure
    eos = eos
    p_bc = 9.28508e-03
#    value = 1.3720000000000002e-006
    boundary = 'left'
  [../]

  [./RadiationRight]
#    type = DirichletBC
    type = RheaBCs
    variable = epsilon
    equation_name = RADIATION
    velocity = velocity
    density = rho
    pressure = pressure
    eos = eos
    p_bc = 9.28508e-03
#    value = 1.6642117992569650e-006
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
    solve_type = PJFNK
#    petsc_options = '-snes_mf_operator -snes_ksp_ew'
    petsc_options_iname = '-mat_fd_coloring_err'
    petsc_options_value = '1.e-12'
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
  #num_steps = 1000
  end_time = 12
  dt = 1.e-1
  dtmin = 1e-9
  l_tol = 1e-8
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-6
  l_max_its = 50
  nl_max_its = 50
  [./TimeStepper]
    type = FunctionDT
    time_t =  '0.     1.  '
    time_dt = '1e-2   1e-2'
  [../]
[]

##############################################################################################
#                                        OUTPUT                                              #
##############################################################################################
# Define the functions computing the inflow and outflow boundary conditions.                 #
##############################################################################################

[Outputs]
    output_initial = true
    interval = 10
    console = true
    exodus = true
#    postprocessor_screen = false
#    perf_log = true
[]

##############################################################################################
#                                        DEBUG                                               #
##############################################################################################
# Debug                 #
##############################################################################################

#[Debug]
#  show_var_residual_norms = true
#[]
