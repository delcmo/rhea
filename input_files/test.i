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
Ce = 0.

###### Constans #######
speed_of_light = 0. # 300
a = 1.372e-2
sigma_a0 = 1000.
sigma_t0 = 1000.

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
  [../]
[]

[Functions]
  [./rho_fn]
    type = ParsedFunction
    value = sin(x-t)+2
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

  [./MassHyperbolic]
    type = RheaMass
    variable = rho
    rhou = rhou
  [../]

#  [./MassVisc]
#    type = RheaArtificialVisc
#    variable = rho
#    equation_name = CONTINUITY
#    density = rho
#    internal_energy = internal_energy
#    velocity = velocity
#    radiation = epsilon
#  [../]

  [./MassForcing]
    type = RheaForcingTerms
    variable = rho
    equation_name = CONTINUITY
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
        [./InitialCondition]
            type = ConstantIC
            value = 2.
        [../]
   [../]

   [./internal_energy]
      family = LAGRANGE
   [../]

  [./pressure]
    family = LAGRANGE
  [../]

  [./epsilon]
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
    velocity_PPS_name = AverageVelocity
    eos = eos
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

  [./L2ErrorDensity]
    type = ElementL2Error
    variable = rho
    function = rho_fn
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
  end_time = 3.14
  dt = 6.25e-3
  dtmin = 1e-4
  l_tol = 1e-8
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-6
  l_max_its = 50
  nl_max_its = 50
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
  postprocessor_screen = true
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
