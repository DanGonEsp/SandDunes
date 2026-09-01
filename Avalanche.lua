------------------------------------------------------------------------------------------
-- Navier-Stokes equation, 3d
-- Discretization: Vertex-centered, stabilized
------------------------------------------------------------------------------------------
print("Simulation Begin")
-- Load utility scripts (e.g. from from ugcore/scripts)
ug_load_script ("ug_util.lua")
ug_load_script ("util/load_balancing_util.lua")

ug_load_script("util/domain_disc_util.lua")
ug_load_script("util/conv_rates_kinetic.lua")

RequiredPlugins({"Limex", "NavierStokes"})

local myProblem=require("SandDunesConfig")

------------------------------------------------------------------------------------------
-- Split communicator
------------------------------------------------------------------------------------------
local numProc         = util.GetParamNumber("-numProc", 1, "Number of temporal processes")
local simCase	= util.GetParamNumber("-simCase", 1, "Simulation Case in Table in")-1
local simCaseBnd	= util.GetParamNumber("-simCaseBnd", 1, "Simulation Case (Boundary Conditions)")

SpaceTimeComm = SpaceTimeCommunicator()
SpaceTimeComm:split(numProc)

local rank = ProcRank()
local rank_t=SpaceTimeComm:get_temporal_rank()

local TemporalSize = SpaceTimeComm:get_temporal_size()
local SpaceSize = SpaceTimeComm:get_spatial_size()

print("TemporalSize = " ..TemporalSize)
print("SpaceSize = " ..SpaceSize)

if numProc > 1 then
	simCase = rank_t
elseif simCase < 0 then
	print ("Simulation case (simCase) not available."); exit();
end

------------------------------------------------------------------------------------------
-- Input parameter table
------------------------------------------------------------------------------------------
local csvfile = require "simplecsv"
local InValues, num_rows, num_cols = csvfile.read('./AvalancheTable_in.csv') -- read file csv1.txt to matrix m
if( TemporalSize > num_rows-1) then print ("TemporalSize larger than rows in input parametrs."); exit(); end
if( simCase+1 > num_rows-1) then print ("Simulation case larger than rows in input parametrs."); exit(); end


inflow       = InValues[simCase+2][1]
H_0          = InValues[simCase+2][2]
W0           = InValues[simCase+2][3]
SlipVelValue = InValues[simCase+2][4]

print("Inflow = " ..inflow.."m/s")
print("Heigh = " ..H_0.. "m.")
print("Width = " ..W0.. "m.")
print("SlipVel = " ..SlipVelValue.. "m.")

local fixedNum = string.format("%04d", simCase+1)


------------------------------------------------------------------------------------------
-- parameters
------------------------------------------------------------------------------------------
params =
{
			-- Numerical parameters of the discretization
	dim      = util.GetParamNumber("-dim", 2, "dimensionality of the problem"),
	dir_name = util.GetParam("-dir_name", ""),
	file_name = util.GetParam("-file_name", "Solution"),
	folder_name = util.GetParam("-folder_name", "Solution") .."_".. fixedNum .."_".. "Avanche".. simCaseBnd,
	elem_type = util.GetParam("-elem_type", "quad", "tri, quad"),
	numRefs     = util.GetParamNumber("-numRefs", 4, "number of grid refinements"),
	numPreRefs     = util.GetParamNumber("-numPreRefs", 3, "number of prerefinements (parallel)"),
	
	simCase = simCase,
	simCaseBnd = simCaseBnd,
	
	
	--Output Data
	boolData = util.GetParamBool("-boolData", false),
	data_name = util.GetParam("-data_name", "Data"),
	outputFactor     = util.GetParam("-output", 1, "output every ... steps"),
	boolCheckPoint = util.GetParamBool("-boolCheckPoint", true),
	
	timeMethod = util.GetParam("-timeMethod","limex","euler limex"),
	modifyDT     = util.GetParamBool("-modifyDT", false),
	DT= util.GetParamNumber("-DT", 10.0, "DT[seconds]"),
	DTmin= util.GetParamNumber("-DTmin", 1e-04, "min  DT"),
	numTimeSteps    = util.GetParamNumber("-numTimeSteps", 100, "time steps"),
	
	
	tol     = util.GetParamNumber("-limex-tol", 1e-2, "time step size"),
	nstages = util.GetParamNumber("-limex-nstages", 2, "limex stages (2 default)"),
	limex_partial_mask = util.GetParamNumber("-limex-partial", 0, "limex partial (0 or 3)"),
	limex_debug_level = util.GetParamNumber("-limex-debug-level", 5, "limex debug level (integer)"),
	VelErrorNorm = util.GetParam("-VelErrorNorm","L2","Norm for Pressure error type H1 , L2"),
	PressErrorNorm = util.GetParam("-limexNorm","H1","Norm for Pressure error type H1 , L2"),
	VolErrorNorm = util.GetParam("-VolErrorNorm","L2","Norm for Pressure error type H1 , L2"),
	alphaVel  = util.GetParamNumber("-alphaVel", 1.0, "Error estimator scale factor for Velocity"),
	alphaPress = util.GetParamNumber("-alphaPress", 0.5e-10, "Error estimator scale factor for Pressure"),
	alphaVol = util.GetParamNumber("-alphaVol", 100, "Error estimator scale factor for Volume fraction"),
	
	incr_factor     = util.GetParamNumber("-incr_factor", 1.5),
	red_factor_fail     = util.GetParamNumber("-red_factor_fail", 0.5),
	red_factor_success     = util.GetParamNumber("-red_factor_success", 0.8),
	optimal_newton_steps = util.GetParamNumber("-optimal_newton_steps", 10),
	maxConvRate = util.GetParamNumber("-maxConvRate", 0.9),
	minConvRate = util.GetParamNumber("-minConvRate", 0.5),
	
	max_newton_steps_steady_state=util.GetParamNumber("-max_newton_steps_steady_state", 100),
	max_newton_steps_transient=util.GetParamNumber("-max_newton_steps_transient", 700),
	SteadyAbsDefect = util.GetParamNumber("-AbsDefect", 1e-010),
	SteadyRedDefect = util.GetParamNumber("-RedDefect", 1e-08),
	AbsDefect = util.GetParamNumber("-AbsDefect", 1e-05),
	RedDefect = util.GetParamNumber("-RedDefect", 1e-05),
	NewtonDebug = util.GetParamBool("-NewtonDebug", false),
	NewtonSteadyDebug = util.GetParamBool("-NewtonSteadyDebug", false),
	NewtonUpdater = util.GetParamBool("-NewtonUpdater", true),
	StepDebug = util.GetParamBool("-StepDebug", false),
	
	lambdamaxSteps = util.GetParamNumber("-lambdamaxSteps", 5),
	lambdaStart  = util.GetParamNumber("-lambdaStart", 1.0),

	damping_mg = util.GetParamNumber("-damping_mg", 1.0),
	value_beta = util.GetParamNumber("-value_beta", 0 ),
	--value_beta = util.GetParamNumber("-value_beta", -0.14 ),
	LinAbsDefectImp = util.GetParamNumber("-LinAbsDefectImp", 1e-012),
	LinRedDefectImp = util.GetParamNumber("-LinRedDefectImp", 1e-03),
	LinAbsDefectLim = util.GetParamNumber("-LinAbsDefectLim", 1e-8),
	LinRedDefectLim = util.GetParamNumber("-LinRedDefectLim", 1e-8),
	max_linear_steps_Lim=util.GetParamNumber("-max_linear_steps_lim", 256),
	max_linear_steps_Imp=util.GetParamNumber("-max_linear_steps_imp", 1000),

	
			-- Physical phenomenon of simulation
	doSteadyState = util.GetParamBool("-doSteadyState", false),
	boolSource = util.GetParamBool("-boolSource", false),
	consistentRho_in_source = util.GetParamBool("-consistentRho_in_source", true),
	boolRelativeVel = util.GetParamBool("-boolRelativeVel", true),
	boolGradientPsSource = util.GetParamBool("-boolGradientPsSource", false),
	boolViscPs = util.GetParamBool("-boolViscPs", true),
	boolAveDiff = util.GetParamBool("-boolAveDiff", true),
	boolSlipDiff = util.GetParamBool("-boolSlipDiff", true),
	boolSlipVel = util.GetParamBool("-boolSlipVel", false),
	boolpress_jump= util.GetParamBool("-boolpress_jump", false),
	boolNormal = util.GetParamBool("-boolNormal", false),
	boolFixVel = util.GetParamBool("-boolFixVel", true),
	boolFixVol = util.GetParamBool("-boolFixVol", false),
	boolMassTerm = util.GetParamBool("-boolMassTerm", true),
	boolDensityMean = util.GetParamBool("-boolDensityMean", false),
	
	inflow   = inflow,
	SlipVelValue = SlipVelValue,
	H_0= H_0,
	ReferencePressure  = util.GetParamNumber("-ReferencePressure",  1.7493e2, "interface value"),
	bStokes     = util.GetParamBool("-Stokes", false ,"If defined, only Stokes Eq. computed"),
	bNoLaplace     = util.GetParamNumber("-noLaplace", false,"If defined, only laplace term used"),
	bExactJac     = util.GetParamNumber("-exactJac", 0.0,"If defined, exact jacobian used"),
	bPecletBlend = util.GetParamBool("-PecletBlend", false,"If defined, Peclet Blend used"),
	upwind_m      = util.GetParam("-upwind_m", "full", "Upwind type full or lps"),
	upwind_t      = util.GetParam("-upwind_t", "full", "Upwind type full or lps"),
	upwind_r      = util.GetParam("-upwind_r", "full", "Upwind type full or lps"),
	bPac        = util.GetParamNumber("-pac", false,"If defined, pac upwind used"),
	diffLength  = util.GetParam("-difflength", "raw", "fivepoint, raw, cor Diffusion length type"),
	stab        = util.GetParam("-stab", "fields_2", "Stabilization type (fields or flow viscosity or karimian)"),
	div_correction = util.GetParamBool("-DivCorrection", false ,"Divergence correction for Newton's inner steps'"),
	boolIPVelocity = util.GetParamBool("-boolIPVelocity", true),
	boolTransportJac = util.GetParamBool("-boolTransportJac", true),
	turbViscMethod = util.GetParam("-turbViscMethod","no","TurbVismodel type no , dyn or sma"),
	modellconstant = util.GetParamNumber("-c",0.5),
	update_turb = util.GetParamNumber("-update_turb", 5, "Update Turbulent Viscosity every .. ... iterations"),

	--Material Properties
	nu_a     = util.GetParamNumber("-visc_a", 1.48e-02, "kinematic viscosity"),
	rho_a     = util.GetParamNumber("-rho_a", 1.2, "Air Density"),
	rho_s     = util.GetParamNumber("-rho_s", 2500, "Sand Density"),
	dp     = util.GetParamNumber("-diameter", 1e-03, "Particle Diameter"),
	nu_s     = util.GetParamNumber("-visc_s", 7.104e-09, "kinematic viscosity"),
	c_init        = util.GetParamNumber("-c_init", 1.0, "max volume fraction"),

	alpha_max        = util.GetParamNumber("-alpha_max", 0.635, "max volume fraction"),
	alpha_min        = util.GetParamNumber("-min alpha_min", 0.57, "max volume fraction"),
	packing_factor   = util.GetParamNumber("-packing_factor", 0.6, "Packingfactor"),
	grad_limit = util.GetParamNumber("-grad_limit", 5.0, "grad limit in Normal vector"),
	slope_limit = util.GetParamNumber("-slope_limit", 2.0e-02, "regularization factor in slip and diff velocity"),
	granular_model= util.GetParamNumber("-granular_model", 3, "Opt: 0 Const, 1 Linear, 2 Einstein, 3 Rheology(I) + Einstein, 4 Relax"),
	density_model  = util.GetParam("-density_model", "linear", "constant, linear"),
	drag_mod = util.GetParamNumber("-drag_model", 2, "Opt: 0 StokesLaw, 1 formula, 2 Schiller-Naumann, 3 Turton and Levenspiel"),
	riemman = util.GetParamNumber("-riemman", 2, "Opt: 0 Upwind, 1 Godunov, 2 Rusanov, 3 Roe"),
	--Model 0 pow(0.63+4.8/sqrt(RE),2.0);

	FR = 0.05,
	B_phi = 1,
	deltaGamma = 1e-03,
	Visc_limit = 1e15,

	deltaPs = 1.48e-04,
	deltaI = 1e-03,
	FricMu_1=0.38,
	FricMu_2=0.64,
	I_0 = 0.279,
	gravity = -9.81,
	
}

params.startTime  = 0.0
params.endTime    = params.DT * params.numTimeSteps
params.DTmax = params.DT
params.DTLimex = params.DT

c_init = params.c_init
params.interface_value  = params.alpha_min/params.packing_factor

------------------------------------------------------------------------------------------
-- GridName
------------------------------------------------------------------------------------------

-- Geometry parameters

if not(params.elem_type == "tri" or params.elem_type == "quad") then
	print ("---------------------------------------------------------------------------------------------  ERROR");
	print ("---------------------------------------------------------------------------------------------  ERROR");
	print ("Geometry not found for elemen type = " ..  params.elem_type); exit();
end

params.gridName	= util.GetParam ("-geom","Avalanche"..params.dim.."D_"..params.elem_type..".ugx")


------------------------------------------------------------------------------------------
-- Domain Subsets
------------------------------------------------------------------------------------------

-- Subsets used in the problem
allSubsets = nil
walls = nil
if params.dim == 2 then
	allSubsets = "Inner, Left, Right,Top, Bottom"
	walls = "Left, Right,Top, Bottom"
elseif params.dim == 3 then
	allSubsets = "Inner, Left1, Left2, Right,Top, Bottom, Front1,Front2, Back1, Back2, Back3"
	walls = "Left1, Left2, Right,Top, Bottom, Front1,Front2, Back1, Back2, Back3"
end
Inner_total={"Inner"}



--------------------------------------------------------------------------------
-- Problem setup.
--------------------------------------------------------------------------------

myProblem:Init(params)

------------------------------------------------------------------------------------------
-- FILE NAMES
------------------------------------------------------------------------------------------

SynchronizeProcesses()
vtk_file_name,folder_vtk,folder_name = myProblem:FileNames(rank,SpaceSize)
SynchronizeProcesses()


------------------------------------------------------------------------------------------
-- Initialize UG4
------------------------------------------------------------------------------------------

InitUG (params.dim, AlgebraType("CPU", params.dim+2))



------------------------------------------------------------------------------------------
-- LOG File
------------------------------------------------------------------------------------------

myProblem:LogFiles(rank_t,folder_vtk .. "/LogFile")


------------------------------------------------------------------------------------------
-- Printing Values
------------------------------------------------------------------------------------------


myProblem:PrintingSettings()


------------------------------------------------------------------------------------------
-- load, refine and distribute the grid  (Approximation Space)
------------------------------------------------------------------------------------------

	approxSpace,u = myProblem:ApproximationSpace(allSubsets)

------------------------------------------------------------------------------------------
-- Lua Functions
------------------------------------------------------------------------------------------


---------------------------------------------------------------------- Initial Velocity
--[[function StartValueX(x,y)
	hh=14.1856
	nn=2.5
	cc=math.pow(y/hh,5)
	return params.inflow*(math.min(1.0, math.pow(y/hh,1/nn))*(1-cc) +(cc)* (2*hh - y) * (y ) / (hh * hh))
end]]
function LOGPROF(psi)
	hh=14.1856
	return params.inflow* (2*hh - psi) * (psi ) / (hh * hh)
end

function StartValueX3d(x,y,z)
	return 0.0
end
function StartValueY3d(x,y,z)
	return 0.0
end
function StartValueZ3d(x,y,z)
	return 0.0
end

function StartValueX2d(x,y)
	return 0.0
end
function StartValueY2d(x,y)
	return 0.0
end


---------------------------------------------------------------------- Initial Pressure
Pstd=0.0
function StartValueP2d(x,y)
	return  0.0
end
function StartValueP3d(x,y,z)
	return  0.0
end

---------------------------------------------------------------------- Initial VolumeFraction



function VolumeFraction_1_2d(x,y)
	local a= 0.2
	local b = 0.4
	local c = 0.8
	
	local q = 1.5e-03
	if (x>a and x<b) or (x>c) then
		return q
	else
		return 0.0
	end
end
function VolumeFraction_2_2d(x,y)
	local a= 0.2
	local b = 0.4
	local c = 0.8
	
	local q = 1.5e-03
	if (x>a and x<b) or (x>c) then
		return q
	else
		return 0.0
	end
end
function VolumeFraction_1_3d(x,y,z)
	return 0.0
end


---------------------------------------------------------------------- Boundary Condition
----------------------------------------------------------- Inlet

function InflowVel2d(x, y, t)
	return StartValueX2d(x,y),StartValueY2d(x,y)
end
function MassInflowVel2d(x, y, t)
	local rho = 1.0--params.rho_a
	return  rho * StartValueX2d(x,y), rho * StartValueY2d(x,y)
end

function InflowVel3d(x, y, z, t)
	return StartValueX3d(x,y,z),StartValueY3d(x,y,z),StartValueZ3d(x,y,z)
end
function MassInflowVel3d(x, y, z, t)
	local rho = 1.0--params.rho_a
	return  rho * StartValueX3d(x,y,z), rho * StartValueY3d(x,y,z), rho * StartValueZ3d(x,y,z)
end


----------------------------------------------------------- Top

function TopFlux12d(x,y)
	local a= 0.3
	local b = 2
	local q = -0.001
	if x>a and x<b then
		return q
	else
		return 0.0
	end
end
function TopFlux22d(x,y)
	local a= 0.2
	local b = 0.4
	local c = 0.8
	
	local q = -0.01
	if (x>a and x<b) or (x>c) then
		return q
	else
		return 0.0
	end
end

function TopFlux13d(x,y,z)
	local a= -1
	local b = 2
	local q = -0.001
	if x>a and x<b then
		return q
	else
		return 0.0
	end
end
function TopFlux23d(x,y,z)
	local a= -1
	local b = 2
	local q = -0.001
	if x>a and x<b then
		return q
	else
		return 0.0
	end
end


-------------------------------------------------------------------------- Parameters List
------------------------------------------------------------------------------------------
-- Parameters List
------------------------------------------------------------------------------------------


InterfaceValues = myProblem:InterfaceParameters()


------------------------------------------------------------------------------------------
-- Secondary Variables (Closures)
------------------------------------------------------------------------------------------


myProblem:Clousures(approxSpace,u,walls)


------------------------------------------------------------------------------------------
-- Compose the discretization
------------------------------------------------------------------------------------------


NavierStokesDisc = myProblem:Discretization(Inner_total)


------------------------------------------------------------------------------------------
-- Boundary Conditions
------------------------------------------------------------------------------------------


InletDisc = NavierStokesInflowFV1M (NavierStokesDisc)
if params.dim == 2 then
	InletDisc:add ("InflowVel"..params.dim.."d", "InflowVel"..params.dim.."d","Left,Top")
elseif params.dim == 3 then
	InletDisc:add ("InflowVel"..params.dim.."d", "InflowVel"..params.dim.."d",walls)
	--InletDisc:add ("InflowVel3d", "InflowVel3d","Back,Front")
end



-- boundary condition at the impermeable walls
WallDisc = NavierStokesWall (NavierStokesDisc)
WallDisc:add ("Bottom")

local DirichletBnd = DirichletBoundary()
local NeumannBnd = NeumannBoundaryFV1("c")
NeumannBnd:add("TopFlux"..simCaseBnd..params.dim.."d","Top", "Inner")
DirichletBnd:add(1.0, "c", "Bottom")
if params.dim == 2 then
	if simCaseBnd == 1 then
		DirichletBnd:add(0.0, "c", "Left")
		NeumannBnd:add(0.0,"Right", "Inner")

	elseif simCaseBnd == 2 then
		DirichletBnd:add(0.0, "c", "Left")
		DirichletBnd:add(0.0, "c", "Right")
		
		
	else
		print ("simCaseBnd Not defined"); exit();
	end
elseif params.dim == 3 then

	DirichletBnd:add(1.0, "c", "Bottom")
	if simCaseBnd == 1 then
		DirichletBnd:add(0.0, "c", "Left1,Left2,Front1,Front2,Right,Back1,Back2,Back3")

	else
		DirichletBnd:add(0.0, "c", "Left2,Front1,Front2,Right,Back1,Back2,Back3")
		NeumannBnd:add(0.0,"Left1,Back1", "Inner")
	end

end


--flowBnd:add(0.0, "v", "Right")
--flowBnd:add(0.0, "p", "Right")
--flowBnd:add(ConstValue, "c", "Top")
--flowBnd:add(ConstValue, "c", "Left")
--flowBnd:add(1.0, "c", "Bottom")


---------------------------------------------------------------------------------------
-- Parameters Inputs
---------------------------------------------------------------------------------------
if not(params.bStokes) then
	myProblem.Density:set_volume_fraction(NavierStokesDisc:volume_fraction())
end

myProblem.Diffusion:set_velocity_gradient(NavierStokesDisc:velocity_grad())

--myProblem.DensityRelVel:set_viscosity(NavierStokesDisc:einstein_viscosity())

myProblem.KinMixViscosity:set_import_2(NavierStokesDisc:mix_viscosity())

---------------------------------------------------------------------------------------
-- Global Discretization
---------------------------------------------------------------------------------------

-- the global discretization
domainDisc = DomainDiscretization (approxSpace)
domainDisc:add (NavierStokesDisc)
domainDisc:add (InletDisc)
domainDisc:add (WallDisc)
domainDisc:add(DirichletBnd)
domainDisc:add(NeumannBnd)



--domainDisc:add(TransportEq)
--domainDisc:add(OutflowBND)

print("Domain Discretization: DONE")

---------------------------------------------------------------------------------------
-- Time Discretization
---------------------------------------------------------------------------------------

print("Setting Time Discretization")

local timeDisc = myProblem:TimeDiscretization(domainDisc)


------------------------------------------------------------------------------------------
-- Interpolate initial values
------------------------------------------------------------------------------------------
print("Initializing Values")
-- start
time = 0
step = 0
local time_work_total = 0.0
local interpolate = false

if(params.boolCheckPoint) then
	time, step, time_work_total, interpolate = myProblem:LoadCheckPoint(u,folder_vtk)
end

if interpolate then
	--Interpolate(StartValueX, u, "u")
	Interpolate(0.0, u, "u")
	Interpolate("StartValueY"..params.dim.."d", u, "v")
	if params.dim == 3 then
		Interpolate("StartValueZ"..params.dim.."d", u, "v")
	end
	Interpolate("StartValueP"..params.dim.."d", u, "p")
	Interpolate("VolumeFraction_"..simCaseBnd.."_"..params.dim.."d", u, "c")
	print("Initial Conditions: Done")
end

------------------------------------------------------------------------------------------
-- Set up the solver
------------------------------------------------------------------------------------------
print("Setting Solver")
boolSolution = 1
op, NLSolver, NewtonSolverSteady, limex = myProblem:CreateSolver(domainDisc, approxSpace)



------------------------------------------------------------------------------------------
-- Set up the Output (For printing variables)
------------------------------------------------------------------------------------------


out = myProblem:OutputParameters()




myProblem.KinTurbulentViscosity:update()
myProblem.gamma:update()
myProblem.RelVel:update()
if params.boolSlipDiff then
	myProblem.SlipDiff:update()
else if params.boolSlipVel then
		myProblem.SlipVel:update()
	end
end
myProblem.Normal:update()

------------------------------------------------------------------------------------------
-- Steady State Solution
------------------------------------------------------------------------------------------
print("Calculating SteadyState")
local time_work_steady=0.0
local linsolver_calls = 0
local linsolver_steps = 0

if params.doSteadyState and interpolate then
	-- Steady state solution.
	
	NewtonSolverSteady:add_inner_step_update(myProblem.gamma)
	NewtonSolverSteady:add_step_update(myProblem.RelVel)
	NewtonSolverSteady:add_step_update(myProblem.Normal)
	if params.turbViscMethod=="no" then
		NewtonSolverSteady:add_step_update(myProblem.KinTurbulentViscosity)
	else
		NewtonSolverSteady:add_inner_step_update(myProblem.KinTurbulentViscosity)
	end

	if params.boolSlipDiff then
		NewtonSolverSteady:add_step_update(myProblem.SlipDiff)
	else if params.boolSlipVel then
			NewtonSolverSteady:add_step_update(myProblem.SlipVel)
		end
	end
	
	time_work_steady, linsolver_calls, linsolver_steps, boolSolution = myProblem:ComputeNonLinearSteadyStateSolution(u, domainDisc, NewtonSolverSteady)
	time_work_total = time_work_total + time_work_steady
end

if(params.boolFixVel) then
	fixer = DirichletBoundary()
	domainDisc:add(fixer)
	fixer:invert_subset_selection()
	fixer:add("u", "")
	fixer:add("v", "")
	fixer:add("p", "")
end
if(params.boolFixVol) then
	fixer = DirichletBoundary()
	domainDisc:add(fixer)
	fixer:invert_subset_selection()
	fixer:add("c", "")
end
------------------------------------------------------------------------------------------
-- Updating attachments
------------------------------------------------------------------------------------------

if params.turbViscMethod=="no" then
	NLSolver:add_step_update(myProblem.KinTurbulentViscosity)
else
	NLSolver:add_inner_step_update(myProblem.KinTurbulentViscosity)
end


NLSolver:add_step_update(myProblem.RelVel)
NLSolver:add_step_update(myProblem.Normal)

if params.timeMethod == "limex" then
	NLSolver:add_step_update(myProblem.gamma)
	if params.boolSlipDiff then
		NLSolver:add_step_update(myProblem.SlipDiff)
	elseif params.boolSlipVel then
		NLSolver:add_step_update(myProblem.SlipVel)
	end
else
	NLSolver:add_inner_step_update(myProblem.gamma)
	NLSolver:add_inner_step_update(myProblem.Normal)
	if params.boolSlipDiff then
		NLSolver:add_inner_step_update(myProblem.SlipDiff)
	elseif params.boolSlipVel then
		NLSolver:add_inner_step_update(myProblem.SlipVel)
	end
	
end
	

------------------------------------------------------------------------------------------
-- Printing Initial Conditions
------------------------------------------------------------------------------------------

	-- write start solution
if boolSolution == 1 then

	print("Writing initial values")
	out:print_subsets(vtk_file_name, u,allSubsets,step,time, true)
	print ("Output to file " .. vtk_file_name .. ".vtu  in time t = 0")
	print ("    -   -   -   -   -   -   -   -   -   -   -   -   -   -   ")
	print ("                                                            ")
	print ("    -   -   -   -   -   -   -   -   -   -   -   -   -   -   ")
	
	myProblem:SaveCheckPoint(u,folder_vtk)
	
	Value_inner1 = Integral(NavierStokesDisc:volume_fraction(), u,"Inner",0.0)

	if (rank == 0 and interpolate) then
		myProblem:WriteValues( folder_vtk, step, time, Value_inner1, 0.0, time_work_steady, time_work_total, 1, 0, linsolver_calls, linsolver_steps,false)
	end
	
end


------------------------------------------------------------------------------------------
-- Final Setting
------------------------------------------------------------------------------------------
-- create new grid function for old value
uOld = u:clone()

-- store grid function in vector of  old solutions
solTimeSeries = SolutionTimeSeries()
solTimeSeries:push(uOld, time)






total_Newton_Steps = 0
total_Newton_Steps_fail = 0
total_linsolver_calls_step = 0
total_linsolver_steps_step = 0

--doo = true
------------------------------------------------------------------------------------------
-- Time Steps Loop    (Solution)
------------------------------------------------------------------------------------------
if boolSolution == 1 then
	for step = step+1, params.numTimeSteps do

		print("++++++ TIMESTEP " .. step .. " BEGIN ++++++")
		tBefore_step = os.clock()
		StartTime = time
		EndTime = time + params.DT
		if params.timeMethod == "limex" then
		
			Newton_Steps, Newton_Steps_fail, linsolver_calls_step, linsolver_steps_step, boolSolution  = myProblem:SolveNonlinearProblemLimex(u, limex, NLSolver, step, StartTime, EndTime, NewtonLimexSteps)

		else
			--[[if doo then
				
				for step2 = 1, 1 do
					Newton_Steps2, Newton_Steps_fail2, linsolver_calls_step2, linsolver_steps_step2 , boolSolution = myProblem:SolveNonlinearProblem( u, NLSolver, op, solTimeSeries, 1, 0,0,1)
				
				end
				fixer = DirichletBoundary()
				domainDisc:add(fixer)
				fixer:invert_subset_selection()
				fixer:add("c", "")
				doo = false
			end]]
			Newton_Steps, Newton_Steps_fail, linsolver_calls_step, linsolver_steps_step , boolSolution = myProblem:SolveNonlinearProblem( u, NLSolver, op, solTimeSeries, DT, step,StartTime,EndTime)
		end
		time = EndTime
		tAfter_step = os.clock()
		
		
		if boolSolution == 1 then
		
			if (step % params.outputFactor == 0 ) then
				out:print_subsets(vtk_file_name, u,allSubsets,step,time)
				print ("Output to file " .. vtk_file_name .. ".vtu  in time t =  " .. time .. "  Step = " .. step .. ".")
				print(" ")
			end
			
			myProblem:SaveCheckPoint(u,folder_vtk)
			
			print("++++++ TIMESTEP " .. step .. "  END ++++++")
			print ("    -   -   -   -   -   -   -   -   -   -   -   -   -   -   ")
			print ("                                                            ")
			print ("                                                            ")
			print ("                                                            ")
			print ("    -   -   -   -   -   -   -   -   -   -   -   -   -   -   ")
			print ("                                                            ")
			print ("<<<<<< Total Newton semi   Steps =  " .. Newton_Steps .. "   >>>>>>")
			print ("<<<<<<       Newton success Steps =  " .. Newton_Steps-Newton_Steps_fail .. "   >>>>>>")
			print ("<<<<<<       Newton fail   Steps =  " .. Newton_Steps_fail .. "   >>>>>>")
			print ("                                                            ")
			print ("    -   -   -   -   -   -   -   -   -   -   -   -   -   -   ")
			print ("                                                            ")
			print ("                                                            ")
			print ("                                                            ")
			print ("    -   -   -   -   -   -   -   -   -   -   -   -   -   -   ")
			
			total_Newton_Steps = total_Newton_Steps + Newton_Steps
			total_Newton_Steps_fail = total_Newton_Steps_fail + Newton_Steps_fail
			total_linsolver_calls_step = total_linsolver_calls_step + linsolver_calls_step
			total_linsolver_steps_step = total_linsolver_steps_step + linsolver_steps_step
			
					
			Value_inner1 = Integral(NavierStokesDisc:volume_fraction(), u,"Inner",0.0)
			
			if rank == 0 and boolSolution == 1 then
				time_work_total = time_work_total + tAfter_step - tBefore_step
				myProblem:WriteValues( folder_vtk, step, time, Value_inner1, 0.0, tAfter_step - tBefore_step, time_work_total, Newton_Steps, Newton_Steps_fail, linsolver_calls_step, linsolver_steps_step,false)
			end
			
		else
			print("++++++ TIMESTEP " .. step .. "  FAILED ++++++")
			out:print_subsets(vtk_file_name, u,allSubsets,step,time)
			print ("Failed Output file" .. vtk_file_name .. ".vtu  in time t =  " .. time .. "  Step = " .. step .. ".")
			print("++++++ TIMESTEP " .. step .. "  FAILED ++++++")
			print(" ")
			break
		end
		
		
	end
end


if boolSolution == 1 then
	------------------------------------------------------------------------------------------
	-- Solution Done
	------------------------------------------------------------------------------------------

	print("-			-")
	print("-------------------------------------------------------------------------------")
	print("Steady state Computation took " .. time_work_steady .. " seconds.")
	print("Temporal Computation took " .. time_work_total-time_work_steady .. " seconds.")
	print("Total Computation took " .. time_work_steady .. " seconds.")
	print("-------------------------------------------------------------------------------")
	print("")
	print("")
	print ("Output to file " .. vtk_file_name .. ".vtu")
	print("done.")


	if rank == 0 then
		myProblem:WriteValues( folder_vtk, params.numTimeSteps, time, Value_inner1, 0.0, time_work_total, time_work_total, total_Newton_Steps, total_Newton_Steps_fail, total_linsolver_calls_step, total_linsolver_steps_step,true)
	end
end

SynchronizeProcesses()
if (params.NewtonDebug and rank == 0 and SpaceSize > 1) then

	csvfile.MergeDebugPVD(params.debug_dir, params.file_name )
	print("NewtonDebug Done")
	
end

--[[local Tablename = folder_name .. "/Table_out_" .. numProc ..".csv"
lineWriter = LineWriter()
Headers = " Sim, Vel, H0, W0, Solved\n"
lineWriter:write_line(Tablename,simCase, Headers, params.inflow, H_0, W0, boolSolution)
]]

myProblem:RunParaViewContour( rank, folder_vtk)

--SynchronizeProcesses()
--SpaceTimeComm:unsplit()

