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
sigma_a0 = 4.2e+004 # 3.9071164263502113e+004
sigma_t0 = 8.5314410158161813e+000

###### Initial Conditions #######
rho_init_left = 1.
rho_init_right = 3.5979106530611014
vel_init_left = 5.8566201857385288e+001
vel_init_right = 1.6277836640428294e+001 
temp_init_left = 1.
temp_init_right = 8.5571992184757608e+000
p_init_left = 82.32
p_init_right = 390.9772444
eps_init_left = 1.3720000000000000e-002
eps_init_right = 7.3566599630083758e-001
membrane = 0.
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
  nx = 2000
  xmin = -9.040524501451145667e-2
  xmax = 2.e-2
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

    [./jump_grad_dens]
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
    jump_dens = jump_grad_dens
    epsilon_PPS_name = AverageEpsilon
    velocity_PPS_name = AverageVelocity
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
    p_bc = 2534.471307
    boundary = 'right'
  [../]

  [./MomentumLeft]
    type = DirichletBC
    variable = rhou
    value = 5.8566201857385288e+001
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
    p_bc = 2534.471307
    boundary = 'right'
  [../]

  [./EnergyLeft]
    type = DirichletBC
    variable = rhoE
    value = 1.8384800000000000e+003
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
    p_bc = 2534.471307
    boundary = 'right'
  [../]

  [./RadiationLeft]
#    type = DirichletBC
#    value = 1.3720000000000000e-002
    type = RheaBCs
    variable = epsilon
    equation_name = RADIATION
    velocity = velocity
    density = rho
    pressure = pressure
    eos = eos
    p_bc = 2534.471307
    boundary = 'left'
  [../]

  [./RadiationRight]
    type = DirichletBC
    variable = epsilon
    value = 7.3566599630083758e+001 
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
#num_steps = 2000
  end_time = 0.05
  dt = 1.e-6
  dtmin = 1e-9
  l_tol = 1e-8
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-6
  l_max_its = 50
  nl_max_its = 50
  [./TimeStepper]
    type = FunctionDT
    time_t =  '0.     1e-3    1.'
    time_dt = '1e-8   1e-4    1e-4'
  [../]
[]

##############################################################################################
#                                        OUTPUT                                              #
##############################################################################################
# Define the functions computing the inflow and outflow boundary conditions.                 #
##############################################################################################

[Output]
  file_base = Mach5_medium_2000_out
  output_initial = true
  interval = 50
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
