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
speed_of_light = 0.
a = 0.
sigma_a0 = 0.
sigma_t0 = 1.

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
  nx = 1280
  xmin = 0
  xmax = 6.28
  block_id = '0'
[]

[Functions]
  [./rho_fn]
    type = ParsedFunction
    value = sin(x-t)+2
  [../]

  [./rhou_fn]
    type = ParsedFunction
    value = (sin(x-t)+2)*(cos(x-t)+2)
  [../]

  [./rhoE_fn]
    type = ParsedFunction
    value = 0.5*1.4/0.4*(cos(x-t)+2)+0.5*(sin(x-t)+2)*(cos(x-t)+2)*(cos(x-t)+2)
  [../]
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
        type = FunctionIC
        function = rho_fn
	[../]
  [../]

  [./rhou]
    family = LAGRANGE
    scaling = 1e+0
	[./InitialCondition]
        type = FunctionIC
        function = rhou_fn
	[../]
  [../]

  [./rhoE]
    family = LAGRANGE
    scaling = 1e+0
	[./InitialCondition]
        type = FunctionIC
        function = rhoE_fn
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

  [./MassForcingTerms]
    type = RheaForcingFunction
    variable = rho
    equation_name = CONTINUITY
    eos = eos
  [../]

  [./MomForcingTerms]
    type = RheaForcingFunction
    variable = rhou
    equation_name = MOMENTUM
    eos = eos
  [../]

  [./EnergyForcingTerms]
    type = RheaForcingFunction
    variable = rhoE
    equation_name = ENERGY
    eos = eos
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

  [./L2ErrorDensity]
    type = ElementL2Error
    variable = rho
    function = rho_fn
  [../]

  [./L2ErrorMomentum]
    type = ElementL2Error
    variable = rhou
    function = rhou_fn
  [../]

  [./L2ErrorEnergy]
    type = ElementL2Error
    variable = rhoE
    function = rhoE_fn
  [../]
[]
##############################################################################################
#                               BOUNDARY CONDITIONS                                          #
##############################################################################################
# Define the functions computing the inflow and outflow boundary conditions.                 #
##############################################################################################
[BCs]
  [./MassLeft]
    type = FunctionDirichletBC
    variable = rho
    function = rho_fn
    boundary = 'left'
  [../]

  [./MassRight]
    type = FunctionDirichletBC
    variable = rho
    function = rho_fn
    boundary = 'right'
  [../]

  [./MomentumLeft]
    type = FunctionDirichletBC
    variable = rhou
    function = rhou_fn
    boundary = 'left'
  [../]

  [./MomentumRight]
    type = FunctionDirichletBC
    variable = rhou
    function = rhou_fn
    boundary = 'right'
  [../]

  [./EnergyLeft]
    type = FunctionDirichletBC
    variable = rhoE
    function = rhoE_fn
    boundary = 'left'
  [../]

  [./EnergyRight]
    type = FunctionDirichletBC
    variable = rhoE
    function = rhoE_fn
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
  #num_steps = 100
  end_time = 1.2 # 3.14
  dt = 0.0015125
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
#postprocessor_screen = false
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
