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

SpaceTimeComm = SpaceTimeCommunicator()
SpaceTimeComm:split(numProc)

local rank = ProcRank()
local rank_t=SpaceTimeComm:get_temporal_rank()

local TemporalSize = SpaceTimeComm:get_temporal_size()
local SpaceSize = SpaceTimeComm:get_spatial_size()

print("TemporalSize = " ..TemporalSize)
print("SpaceSize = " ..SpaceSize)

------------------------------------------------------------------------------------------
-- Input parameter table
------------------------------------------------------------------------------------------
local csvfile = require "simplecsv"
local InValues, num_rows, num_cols = csvfile.read('./FlowTable_in.csv') -- read file csv1.txt to matrix m
if( TemporalSize > num_rows-1) then print ("TemporalSize larger than rows in input parametrs."); exit(); end
if( simCase+1 > num_rows-1) then print ("Simulation case larger than rows in input parametrs."); exit(); end

if numProc > 1 then
	simCase = rank_t
end

inflow       = InValues[simCase+2][1]
H_0          = InValues[simCase+2][2]
W0           = InValues[simCase+2][3]
SlipVelValue = InValues[simCase+2][4]

print("Inflow = " ..inflow.."m/s")
print("Heigh = " ..H_0.. "m.")
print("Width = " ..W0.. "m.")
print("SlipVel = " ..SlipVelValue.. "m.")

local fixedNum = string.format("%04d", simCase)


------------------------------------------------------------------------------------------
-- parameters
------------------------------------------------------------------------------------------
params =
{
			-- Numerical parameters of the discretization
	dim      = util.GetParamNumber("-dim", 2, "dimensionality of the problem"),
	file_name = util.GetParam("-file_name", "SolutionX") .."_".. fixedNum,
	folder_name = util.GetParam("-folder_name", "SolutionY") .. "Flow",
	elem_type = util.GetParam("-elem_type", "quad", "tri, quad"),
	numRefs     = util.GetParamNumber("-numRefs", 3, "number of grid refinements"),
	numPreRefs     = util.GetParamNumber("-numPreRefs", 1, "number of prerefinements (parallel)"),
	
	simCase = simCase,
	
	DT= util.GetParamNumber("-DT", 1.0, "DT[seconds]"),
	DTmin= util.GetParamNumber("-DTmin", 1e-04, "min  DT"),
	numTimeSteps    = util.GetParamNumber("-numTimeSteps", 10, "time steps"),
	outputFactor     = util.GetParam("-output", 1, "output every ... steps"),
	
	timeMethod = util.GetParam("-timeMethod","euler","euler limex"),
	modifyDT     = util.GetParamBool("-modifyDT", false),
	
	tol     = util.GetParamNumber("-limex-tol", 1e-1, "time step size"),
	nstages = util.GetParamNumber("-limex-nstages", 2, "limex stages (2 default)"),
	limex_partial_mask = util.GetParamNumber("-limex-partial", 0, "limex partial (0 or 3)"),
	limex_debug_level = util.GetParamNumber("-limex-debug-level", 5, "limex debug level (integer)"),
	VelErrorNorm = util.GetParam("-VelErrorNorm","L2","Norm for Pressure error type H1 , L2"),
	PressErrorNorm = util.GetParam("-limexNorm","H1","Norm for Pressure error type H1 , L2"),
	VolErrorNorm = util.GetParam("-VolErrorNorm","L2","Norm for Pressure error type H1 , L2"),
	alphaVel  = util.GetParamNumber("-alphaVel", 1.0, "Error estimator scale factor for Velocity"),
	alphaPress = util.GetParamNumber("-alphaPress", 0.5e-10, "Error estimator scale factor for Pressure"),
	alphaVol = util.GetParamNumber("-alphaVol", 100, "Error estimator scale factor for Volume fraction"),
	
	incr_factor     = util.GetParamNumber("-incr_factor", 2.0),
	red_factor_fail     = util.GetParamNumber("-red_factor_fail", 0.5),
	red_factor_success     = util.GetParamNumber("-red_factor_success", 0.8),
	optimal_newton_steps = util.GetParamNumber("-optimal_newton_steps", 10),
	maxConvRate = util.GetParamNumber("-maxConvRate", 0.9),
	minConvRate = util.GetParamNumber("-minConvRate", 0.5),
	
	max_newton_steps_steady_state=util.GetParamNumber("-max_newton_steps_steady_state", 100),
	max_newton_steps_transient=util.GetParamNumber("-max_newton_steps_transient", 2500),
	SteadyAbsDefect = util.GetParamNumber("-AbsDefect", 1e-010),
	SteadyRedDefect = util.GetParamNumber("-RedDefect", 1e-08),
	AbsDefect = util.GetParamNumber("-AbsDefect", 1e-010),
	RedDefect = util.GetParamNumber("-RedDefect", 1e-05),
	NewtonDebug = util.GetParamBool("-NewtonDebug", false),
	NewtonSteadyDebug = util.GetParamBool("-NewtonSteadyDebug", false),
	NewtonUpdater = util.GetParamBool("-NewtonUpdater", true),
	StepDebug = util.GetParamBool("-StepDebug", false),
	
	lambdamaxSteps = util.GetParamNumber("-lambdamaxSteps", 5),
	lambdaStart  = util.GetParamNumber("-lambdaStart", 1.0),

	max_linear_steps=util.GetParamNumber("-max_linear_steps", 1000),
	damping_mg = util.GetParamNumber("-damping_mg", 0.9),
	value_beta = util.GetParamNumber("-value_beta", -0.10 ),
	--value_beta = util.GetParamNumber("-value_beta", -0.14 ),
	LinAbsDefectImp = util.GetParamNumber("-LinAbsDefectImp", 1e-012),
	LinRedDefectImp = util.GetParamNumber("-LinRedDefectImp", 1e-04),
	LinAbsDefectLim = util.GetParamNumber("-LinAbsDefectLim", 1e-018),
	LinRedDefectLim = util.GetParamNumber("-LinRedDefectLim", 1e-12),

	
			-- Physical phenomenon of simulation
	doSteadyState = util.GetParamBool("-doSteadyState", true),
	boolSource = util.GetParamBool("-boolSource", false),
	consistentRho_in_source = util.GetParamBool("-consistentRho_in_source", true),
	boolRelativeVel = util.GetParamBool("-boolRelativeVel", true),
	boolGradientPsSource = util.GetParamBool("-boolGradientPsSource", false),
	boolViscPs = util.GetParamBool("-boolViscPs", true),
	boolAveDiff = util.GetParamBool("-boolAveDiff", true),
	boolSlipDiff = util.GetParamBool("-boolSlipDiff", false),
	boolSlipVel = util.GetParamBool("-boolSlipVel", true),
	boolpress_jump= util.GetParamBool("-boolpress_jump", false),
	boolNormal = util.GetParamBool("-boolNormal", false),
	boolFixVel = util.GetParamBool("-boolFixVel", false),
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
	lee_factor = util.GetParamNumber("-lee_factor", 1.0, "lee slope factor"),
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

params.gridName	= util.GetParam ("-geom","Dune"..params.dim.."D_"..params.elem_type.."_double.ugx")


------------------------------------------------------------------------------------------
-- Domain Subsets
------------------------------------------------------------------------------------------

-- Subsets used in the problem
allSubsets = "Inner, Inner2, Left, Right,Top, Bottom"
if params.dim == 3 then
	allSubsets = allSubsets .. ", Back, Front"
end
Inner_total={"Inner","Inner2"}



--------------------------------------------------------------------------------
-- Problem setup.
--------------------------------------------------------------------------------
local myProblem=require("SandDunesConfig")
myProblem:Init(params)

------------------------------------------------------------------------------------------
-- FILE NAMES
------------------------------------------------------------------------------------------

SynchronizeProcesses()
vtk_file_name,folder,folder_name = myProblem:FileNames(rank)
SynchronizeProcesses()


------------------------------------------------------------------------------------------
-- Initialize UG4
------------------------------------------------------------------------------------------

InitUG (params.dim, AlgebraType("CPU", params.dim+2))
GetLogAssistant():enable_file_output(true, folder .. "/LogFile_"..simCase.. "_Lev"..params.numRefs)
if rank_t > 0 then GetLogAssistant():enable_terminal_output(false) end


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
	return LOGPROF(z)
end
function StartValueY3d(x,y,z)
	return 0.0
end
function StartValueZ3d(x,y,z)
	return 0.0
end

function StartValueX2d(x,y)
	return LOGPROF(y)
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
function PressureBoundary(x,y)
	return  1.2*y*params.gravity
end

---------------------------------------------------------------------- Initial VolumeFraction

sigma=3
mu_c=20

---------------------------------------------------------------------- Initial DuneShape
-- Super-Gaussian parameters
--H_0 = 20       -- Dune heigh
x0 = 20          -- Center in x
y0 = 0           -- Center in y
sigmaX = 3.0     -- Width in x
sigmaY = 3.0     -- Width in y
n = 4            -- Super-Gaussian order (2 = Gaussian)
h_0 = 0.01         -- offset

-- Super-Gaussian function
function superGaussian3d(x, y)
    local dx = math.abs((x - x0) / sigmaX)
    local dy = math.abs((y - y0) / sigmaY)

    return H_0 * math.exp(-(dx^n + dy^n))- h_0
end
function superGaussian2d(x)
    local dx = math.abs((x - x0) / sigmaX)
    return H_0 * math.exp(-(dx^n)) - h_0
end


function VolumeFraction2d(x,y)
	if y>superGaussian2d(x) then
		return 0.0
	else
		return c_init
	end
end
function VolumeFraction3d(x,y,z)
	if z>superGaussian3d(x,y) then
		return 0.0
	else
		return c_init
	end
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


----------------------------------------------------------- Bottom

function BoundaryVolumeFraction(x,y)
	a= 10
	b = 20
	if x>a and x<b then--y>Dune(x,y) then
		return c_init
	else
		return 0.0
	end
end

function BottomFlux(x,y) return 0 end
-------------------------------------------------------------------------- Parameters List
------------------------------------------------------------------------------------------
-- Parameters List
------------------------------------------------------------------------------------------


InterfaceValues = myProblem:InterfaceParameters()


------------------------------------------------------------------------------------------
-- Secondary Variables (Closures)
------------------------------------------------------------------------------------------


myProblem:Clousures(approxSpace,u)


------------------------------------------------------------------------------------------
-- Compose the discretization
------------------------------------------------------------------------------------------


NavierStokesDisc = myProblem:Discretization(Inner_total)


------------------------------------------------------------------------------------------
-- Boundary Conditions
------------------------------------------------------------------------------------------


InletDisc = NavierStokesInflowFV1M (NavierStokesDisc)
InletDisc:add ("InflowVel"..params.dim.."d", "InflowVel"..params.dim.."d","Left,Top")
if params.dim == 3 then
	InletDisc:add ("InflowVel3d", "InflowVel3d","Back,Front")
end


-- boundary condition at the outlet
OutletDisc = NavierStokesNoNormalStressOutflowFV1M (NavierStokesDisc)
OutletDisc:add ("Right")
OutletDisc:set_phase_parameters(myProblem.InterfaceValues)

-- boundary condition at the impermeable walls
WallDisc = NavierStokesWall (NavierStokesDisc)
WallDisc:add ("Bottom")


flowBnd = DirichletBoundary()
flowBnd:add(0.0, "c", "Left,Top")
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
domainDisc:add (OutletDisc)
domainDisc:add (WallDisc)
domainDisc:add(flowBnd)


--domainDisc:add(TransportEq)
--domainDisc:add(OutflowBND)

print("Domain Discretization: DONE")
---------------------------------------------------------------------------------------
-- Time Discretization
---------------------------------------------------------------------------------------
print("Time Discretization")

local timeDisc = myProblem:TimeDiscretization(domainDisc)

------------------------------------------------------------------------------------------
-- Set up the solver
------------------------------------------------------------------------------------------
print("Setting Solver")
boolSolution = 1
op, NLSolver, NewtonSolverSteady, limex, boolSolution = myProblem:CreateSolver(domainDisc, approxSpace)



------------------------------------------------------------------------------------------
-- Set up the Output (For printing variables)
------------------------------------------------------------------------------------------


out = myProblem:OutputParameters()


------------------------------------------------------------------------------------------
-- Interpolate initial values
------------------------------------------------------------------------------------------
print("Initializing Values")

--Interpolate(StartValueX, u, "u")
Interpolate(0.0, u, "u")
Interpolate("StartValueY"..params.dim.."d", u, "v")
if params.dim == 3 then
	Interpolate("StartValueZ"..params.dim.."d", u, "v")
end
Interpolate("StartValueP"..params.dim.."d", u, "p")
Interpolate("VolumeFraction"..params.dim.."d", u, "c")


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
time_work_steady=0.0
linsolver_calls = 0
linsolver_steps = 0

if params.doSteadyState and boolSolution == 1 then
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
-- Printing Initial Conditions
------------------------------------------------------------------------------------------
-- start
time = 0
step = 0

	-- write start solution
if boolSolution == 1 then

	print("Writing initial values")
	out:print_subsets(vtk_file_name, u,allSubsets,step,time, true)
	print ("Output to file " .. vtk_file_name .. ".vtu  in time t = 0")
	print ("    -   -   -   -   -   -   -   -   -   -   -   -   -   -   ")
	print ("                                                            ")
	print ("    -   -   -   -   -   -   -   -   -   -   -   -   -   -   ")
	
end


------------------------------------------------------------------------------------------
-- Final Setting
------------------------------------------------------------------------------------------
-- create new grid function for old value
uOld = u:clone()

-- store grid function in vector of  old solutions
solTimeSeries = SolutionTimeSeries()
solTimeSeries:push(uOld, time)


Value_inner1 = Integral(NavierStokesDisc:volume_fraction(), u,"Inner",0.0)
Value_inner2 = Integral(NavierStokesDisc:volume_fraction(), u,"Inner2",0.0)

if (rank == 0 and  boolSolution == 1) then
	myProblem:WriteValues( folder, step, time, Value_inner1, Value_inner2, time_work_steady, 1, 0, linsolver_calls, linsolver_steps,false)
end


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
	if params.boolSlipDiff then
		NLSolver:add_inner_step_update(myProblem.SlipDiff)
	elseif params.boolSlipVel then
		NLSolver:add_inner_step_update(myProblem.SlipVel)
	end
end
	



total_Newton_Steps = 0
total_Newton_Steps_fail = 0
total_linsolver_calls_step = 0
total_linsolver_steps_step = 0
tBefore = os.clock()

--doo = true
------------------------------------------------------------------------------------------
-- Time Steps Loop    (Solution)
------------------------------------------------------------------------------------------
if boolSolution == 1 then
	for step = 1, params.numTimeSteps do

		print("++++++ TIMESTEP " .. step .. " BEGIN ++++++")
		tBefore_step = os.clock()
		StartTime = time
		EndTime = time + params.DT
		if params.timeMethod == "limex" then
		
			Newton_Steps, Newton_Steps_fail, linsolver_calls_step, linsolver_steps_step, boolSolution  = myProblem:SolveNonlinearProblemLimex(u, limex, NLSolver, step, StartTime, EndTime)

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
			Value_inner2 = Integral(NavierStokesDisc:volume_fraction(), u,"Inner2",0.0)
			
			if rank == 0 then
				myProblem:WriteValues( folder, step, time, Value_inner1, Value_inner2, tAfter_step - tBefore_step, Newton_Steps, Newton_Steps_fail, linsolver_calls_step, linsolver_steps_step,false)
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
tAfter = os.clock()

if boolSolution == 1 then
	------------------------------------------------------------------------------------------
	-- Solution Done
	------------------------------------------------------------------------------------------

	print("-			-")
	print("-------------------------------------------------------------------------------")
	print("Steady state Computation took " .. time_work_steady .. " seconds.")
	print("Temporal Computation took " .. tAfter-tBefore .. " seconds.")
	print("Total Computation took " .. time_work_steady+tAfter-tBefore .. " seconds.")
	print("-------------------------------------------------------------------------------")
	print("")
	print("")
	print ("Output to file " .. vtk_file_name .. ".vtu")
	print("done.")


	if rank == 0 then
		myProblem:WriteValues( folder, params.numTimeSteps, time, Value_inner1, Value_inner2, time_work_steady+tAfter-tBefore, total_Newton_Steps, total_Newton_Steps_fail, total_linsolver_calls_step, total_linsolver_steps_step,true)
	end
end

SynchronizeProcesses()
if (params.NewtonDebug and rank == 0 and SpaceSize > 1) then

	csvfile.MergeDebugPVD(params.debug_dir, params.file_name )
	print("NewtonDebug Done")
	
end

local Tablename = folder_name .. "/Table_out_" .. numProc ..".csv"
lineWriter = LineWriter()
Headers = " Sim, Vel, H0, W0, Solved\n"
lineWriter:write_line(Tablename,simCase, Headers, params.inflow, H_0, W0, boolSolution)

--SynchronizeProcesses()
--SpaceTimeComm:unsplit()

