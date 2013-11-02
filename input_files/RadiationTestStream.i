#
#####################################################
# Define some global parameters used in the blocks. #
#####################################################
#
[GlobalParams]
###### Other parameters #######
order = FIRST

###### Constans #######
speed_of_light = 0.1
a = 1.
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
[]

[Functions]

  [./epsilon_fn]
    type = ParsedFunction
    value=sin(x-1000*t)+2
  [../]

[]

###### Mesh #######
[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 20
  xmin = 0
  xmax = 6.28318
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
        type = FunctionIC
        function = epsilon_fn
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

  [./RadiationForcing]
    type = RheaForcingFunctionStream
    variable = epsilon
    equation_name = RADIATION
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
            value = 2
        [../]
   [../]

   [./temperature]
    family = LAGRANGE
   [../]

  [./pressure]
    family = LAGRANGE
  [../]

  [./rho]
    family = LAGRANGE
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
    eos = eos
  [../]
[]

##############################################################################################
#                                     PPS                                                    #
##############################################################################################
# Define functions that are used in the kernels and aux. kernels.                            #
##############################################################################################
[Postprocessors]

  [./L2ErrorEpsilon]
    type = ElementL2Error
    variable = epsilon
    function = epsilon_fn
  [../]

  [./L1ErrorEpsilon]
    type = ElementH1SemiError
    variable = epsilon
    function = epsilon_fn
  [../]

[]
##############################################################################################
#                               BOUNDARY CONDITIONS                                          #
##############################################################################################
# Define the functions computing the inflow and outflow boundary conditions.                 #
##############################################################################################
[BCs]

   [./Periodic]
    [./EPsilon]
    variable = epsilon
    auto_direction = 'x'
    [../]
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
  end_time = 1e-3.
  dt = 4.e-4
  dtmin = 1e-8
  l_tol = 1e-8
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-6
  l_max_its = 50
  nl_max_its = 50
  [./Quadrature]
    type = TRAP
  [../]
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
#postprocessor_screen = true
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
