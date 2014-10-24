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
sigma_a0 = 3.9071164263502113e+000
sigma_t0 = 8.5314410158161813e+004

###### Initial Conditions #######
rho_init_left = 1.
rho_init_right = 6.5189217901173153
vel_init_left = 5.8566201857385295e+02
vel_init_right = 8.9840319830453623e+01 
temp_init_left = 1.
temp_init_right = 8.5515528368625027e+01
p_init_left = 82.32
p_init_right = 45890.85148
eps_init_left = 1.3720000000000000e-02
eps_init_right = 7.3372623010289914e+05 
membrane = 0.20
length = 0.20
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
  	Cv = 1.2348000000000001e+2
  	q_prime = 0.
  [../]

  [./JumpGradPress]
    type = JumpGradientInterface
    variable = pressure
    jump_name = jump_grad_press
  [../]
[]

###### Mesh #######
[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 300
  xmin = 0
  xmax = 0.30
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
    scaling = 1e-4
	[./InitialCondition]
        type = InitialConditions
        eos = eos
	[../]
  [../]

  [./rhou]
    family = LAGRANGE
    scaling = 1e-6
	[./InitialCondition]
        type = InitialConditions
        eos = eos
	[../]
  [../]

  [./rhoE]
    family = LAGRANGE
    scaling = 1e-6
	[./InitialCondition]
        type = InitialConditions
        eos = eos
	[../]
  [../]

  [./epsilon]
    family = LAGRANGE
    scaling = 1e-6
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
    epsilon = epsilon
    jump_press = jump_grad_press
    epsilon_PPS_name = AverageEpsilon
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

[./AverageEpsilon]
    type = ElementAverageValue
    variable = epsilon
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
    value = 1.
    boundary = 'left'
  [../]

  [./MassRight]
    type = RheaBCs
    variable = rho
    equation_name = CONTINUITY
    velocity = velocity
    temperature = temperature
    density = rho
    pressure = pressure
    epsilon = epsilon
    eos = eos
    p_bc = 45890.85148
    boundary = 'right'
  [../]

  [./MomentumLeft]
    type = DirichletBC
    variable = rhou
    value = 5.8566201857385295e+02
    boundary = 'left'
  [../]

  [./MomentumRight]
    type = RheaBCs
    variable = rhou
    equation_name = MOMENTUM
    velocity = velocity
    temperature = temperature
    density = rho
    pressure = pressure
    epsilon = epsilon
    eos = eos
    p_bc = 45890.85148
    boundary = 'right'
  [../]

  [./EnergyLeft]
    type = DirichletBC
    variable = rhoE
    value = 1.7162348000000004e+05
    boundary = 'left'
  [../]

  [./EnergyRight]
    type = RheaBCs
    variable = rhoE
    equation_name = ENERGY
    velocity = velocity
    temperature = temperature
    density = rho
    pressure = pressure
    epsilon = epsilon
    eos = eos
    p_bc = 45890.85148
    boundary = 'right'
  [../]

  [./RadiationLeft]
    type = DirichletBC
    variable = epsilon
    value = 1.3720000000000000e-02
    boundary = 'left'
  [../]

  [./RadiationRight]
    type = DirichletBC
#    type = RheaBCs
    variable = epsilon
#    equation_name = RADIATION
#    velocity = velocity
#    temperature = temperature
#    density = rho
#    pressure = pressure
#    epsilon = epsilon
#    eos = eos
#    p_bc = 45890.85148
    value = 7.3372623010289914e+05
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
#    petsc_options_iname = '-mat_fd_coloring_err  -mat_fd_type  -mat_mffd_type'
#    petsc_options_value = '1.e-12       ds             ds'
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
#num_steps = 1000
  end_time = 0.005
  dt = 1.e-6
  dtmin = 1e-22
  l_tol = 1e-8
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-6
  l_max_its = 50
  nl_max_its = 10
[./TimeStepper]
    type = FunctionDT
    time_t =  '0.    1e-5  1e-3  4e-3'
    time_dt = '1e-7  1e-7  1e-5  1e-5'
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
    console = true
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
