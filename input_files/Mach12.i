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
sigma_a0 = 3.9071164263502112e+002
sigma_t0 = 8.5314410158161809e+002

###### Initial Conditions #######
rho_init_left = 1.
rho_init_right = 1.2973213452231311e+000
vel_init_left = 1.4055888445772469e-001
vel_init_right = 1.0834546504247138e-001 
temp_init_left = 1.0000000000000001e-001
temp_init_right = 1.1947515210501813e-001
eps_init_left = 1.3720000000000002e-006
eps_init_right = 2.7955320762182542e-006 
p_init_left = 8.232e-3
p_init_right = 1.2759408e-2
membrane = 0.010
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
[]

###### Mesh #######
[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 500
  xmin = 0
  xmax = 0.020
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
        type = InitialConditions
        eos = eos
	[../]
  [../]

  [./rhou]
    family = LAGRANGE
    scaling = 1e+2
	[./InitialCondition]
        type = InitialConditions
        eos = eos
	[../]
  [../]

  [./rhoE]
    family = LAGRANGE
    scaling = 1e+2
	[./InitialCondition]
        type = InitialConditions
        eos = eos
	[../]
  [../]

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

  [./diffusion_aux]
    family = MONOMIAL
    order = CONSTANT
  [../]

 [./sigma_a_aux]
   family = MONOMIAL
    order = CONSTANT
  [../]

  [./sigma_t_aux]
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

  [./DiffusionAK]
    type = MaterialRealAux
    variable = diffusion_aux
    property = diffusion
  [../]

  [./SigmaAAK]
    type = MaterialRealAux
    variable = sigma_a_aux
    property = sigma_a
  [../]

  [./SigmaTAK]
    type = MaterialRealAux
    variable = sigma_t_aux
    property = sigma_t
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
    jump = jump_grad_press
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
  [./MassLeft]
    type = RheaBCs
    variable = rho
    equation_name = CONTINUITY
    velocity = velocity
    temperature = temperature
    density = rho
    pressure = pressure
    eos = eos
    p_bc = 8.232e-3
    T_bc = 1.0000000000000001e-001
    v_bc = 1.4055888445772469e-001
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
    eos = eos
    p_bc = 1.2759408e-2
    boundary = 'right'
  [../]

  [./MomentumLeft]
    type = RheaBCs
    variable = rhou
    equation_name = MOMENTUM
    velocity = velocity
    temperature = temperature
    density = rho
    pressure = pressure
    eos = eos
    p_bc = 8.232e-3
    T_bc = 1.0000000000000001e-001
    v_bc = 1.4055888445772469e-001
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
    eos = eos
    p_bc = 1.2759408e-2
    boundary = 'right'
  [../]

  [./EnergyLeft]
    type = RheaBCs
    variable = rhoE
    equation_name = ENERGY
    velocity = velocity
    temperature = temperature
    density = rho
    pressure = pressure
    eos = eos
    p_bc = 8.232e-3
    T_bc = 1.0000000000000001e-001
    v_bc = 1.4055888445772469e-001
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
    eos = eos
    p_bc = 1.2759408e-2
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
    p_bc = 1.2759408e-2
#    value = 1.3720000000000002e-006
    boundary = 'left'
  [../]

  [./RadiationRight]
    type = DirichletBC
    variable = epsilon
    value = 2.7955320762182542e-006
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
  #num_steps = 1000
  end_time = 1
  dt = 1.e-4
  dtmin = 1e-9
  l_tol = 1e-8
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-6
  l_max_its = 50
  nl_max_its = 10
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
  interval = 100
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
