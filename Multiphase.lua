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
-- Problem
------------------------------------------------------------------------------------------

local problem = util.GetParam("-problem", "avalanche", "flow or avalanche")
local problemTag
if problem == "flow" then problemTag = "MultiphaseFlow" else problemTag = "Avanche" end
if problem ~= "flow" and problem ~= "avalanche" then error("Specify -problem flow or -problem avalanche") end
local defaults = myProblem:GetCaseDefaults(problem)

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
local inputTable
if problem == "flow" then inputTable = "./FlowTable_in.csv" else inputTable = "./AvalancheTable_in.csv" end
local InValues, num_rows, num_cols = csvfile.read(inputTable) -- read file csv1.txt to matrix m
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
			-- Problem configuration
	problem = problem,
	simCase = simCase,
	simCaseBnd = simCaseBnd,
	
	
			-- Numerical parameters of the discretization
	dim      = util.GetParamNumber("-dim", 2, "dimensionality of the problem"),
	dir_name = util.GetParam("-dir_name", ""),
	file_name = util.GetParam("-file_name", "Solution"),
	folder_name = util.GetParam("-folder_name", "Solution") .."_".. fixedNum .."_".. problemTag.. simCaseBnd,
	elem_type = util.GetParam("-elem_type", "quad", "tri, quad"),
	numRefs = util.GetParamNumber("-numRefs", defaults.numRefs, "number of grid refinements"),
	numPreRefs = util.GetParamNumber("-numPreRefs", defaults.numPreRefs, "number of prerefinements (parallel)"),
	
	--Output Data
	boolData = util.GetParamBool("-boolData", false),
	data_name = util.GetParam("-data_name", "Data"),
	outputFactor     = util.GetParam("-output", 1, "output every ... steps"),
	writeIntegral = util.GetParamBool("-writeIntegral", true),
	boolLoadCheckPoint = util.GetParamBool("-boolLoadCheckPoint", true),
	boolSaveCheckPoint = util.GetParamBool("-boolSaveCheckPoint", true),
	
	timeMethod = util.GetParam("-timeMethod","limex","euler limex"),
	modifyDT     = util.GetParamBool("-modifyDT", false),
	DT= util.GetParamNumber("-DT", 1000.0, "DT[seconds]"),
	DTmin= util.GetParamNumber("-DTmin", 1e-04, "min  DT"),
	numTimeSteps    = util.GetParamNumber("-numTimeSteps", 100, "time steps"),
	
	
	tol     = util.GetParamNumber("-limex-tol", 1e-2, "time step size"),
	nstages = util.GetParamNumber("-limex-nstages", 2, "limex stages (2 default)"),
	limex_partial_mask = util.GetParamNumber("-limex-partial", 0, "limex partial (0 or 3)"),
	limex_debug_level = util.GetParamNumber("-limex-debug-level", 5, "limex debug level (integer)"),
	VelErrorNorm = util.GetParam("-VelErrorNorm","L2","Norm for Pressure error type H1 , L2"),
	PressErrorNorm = util.GetParam("-limexNorm","H1","Norm for Pressure error type H1 , L2"),
	VolErrorNorm = util.GetParam("-VolErrorNorm","L2","Norm for Pressure error type H1 , L2"),
	alphaVel = util.GetParamNumber("-alphaVel", defaults.alphaVel, "Error estimator scale factor for Velocity"),
	alphaPress = util.GetParamNumber("-alphaPress", defaults.alphaPress, "Error estimator scale factor for Pressure"),
	alphaVol = util.GetParamNumber("-alphaVol", 100, "Error estimator scale factor for Volume fraction"),
	
	incr_factor     = util.GetParamNumber("-incr_factor", 1.5),
	red_factor_fail     = util.GetParamNumber("-red_factor_fail", 0.7),
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


	------------------------------------------------------------------------------------- LINEAR SOLVER
	damping_mg = util.GetParamNumber("-damping_mg", 1.0),
	rap = util.GetParamBool("-rap", false),
	value_beta = util.GetParamNumber("-value_beta", 0.0 ),
	--value_beta = util.GetParamNumber("-value_beta", -0.14 ),
	LinAbsDefectImp = util.GetParamNumber("-LinAbsDefectImp", 1e-012),
	LinRedDefectImp = util.GetParamNumber("-LinRedDefectImp", defaults.LinRedDefectImp),
	LinAbsDefectLim = util.GetParamNumber("-LinAbsDefectLim", defaults.LinAbsDefectLim),
	LinRedDefectLim = util.GetParamNumber("-LinRedDefectLim", defaults.LinRedDefectLim),
	max_linear_steps_Lim=util.GetParamNumber("-max_linear_steps_lim", 1000),
	max_linear_steps_Imp=util.GetParamNumber("-max_linear_steps_imp", 1000),
	precondLim = util.GetParam("-precondLim","gmg","ilu,gmg"),
	smoother = util.GetParam("-smoother","ilu","ilu,ilut"),
	pre_smooth   = util.GetParamNumber("-pre_smooth", 3, "PreSmooth steps"),
	post_smooth = util.GetParamNumber("-post_smooth", 3, "PostSmooth steps"),
	eps_ilut = util.GetParamNumber("-eps_ilut", 1e-02),


	
			-- Physical phenomenon of simulation
	doSteadyState = util.GetParamBool("-doSteadyState", false),
	boolSource = util.GetParamBool("-boolSource", false),
	consistentRho_in_source = util.GetParamBool("-consistentRho_in_source", true),
	boolRelativeVel = util.GetParamBool("-boolRelativeVel", true),
	boolGradientPsSource = util.GetParamBool("-boolGradientPsSource", false),
	boolViscPs = util.GetParamBool("-boolViscPs", true),
	boolAveDiff = util.GetParamBool("-boolAveDiff", defaults.boolAveDiff),
	boolSlipDiff = util.GetParamBool("-boolSlipDiff", defaults.boolSlipDiff),
	boolSlipVel = util.GetParamBool("-boolSlipVel", defaults.boolSlipVel),
	boolpress_jump= util.GetParamBool("-boolpress_jump", false),
	boolAveNormal = util.GetParamBool("-boolAveNormal", false),
	boolFixVel = util.GetParamBool("-boolFixVel", defaults.boolFixVel),
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
	nu_a     = util.GetParamNumber("-visc_a", 1.48e-05, "kinematic viscosity"),
	rho_a     = util.GetParamNumber("-rho_a", 1.2, "Air Density"),
	rho_s     = util.GetParamNumber("-rho_s", 2500, "Sand Density"),
	dp     = util.GetParamNumber("-diameter", 1e-03, "Particle Diameter"),
	nu_s     = util.GetParamNumber("-visc_s", 1.48e-05, "kinematic viscosity"),
	c_init        = util.GetParamNumber("-c_init", 1.0, "max volume fraction"),

	alpha_max        = util.GetParamNumber("-alpha_max", 0.635, "max volume fraction"),
	alpha_min        = util.GetParamNumber("-min alpha_min", 0.57, "max volume fraction"),
	packing_factor   = util.GetParamNumber("-packing_factor", 0.6, "Packingfactor"),
	grad_limit = util.GetParamNumber("-grad_limit", 0.05, "grad limit in Normal vector"),
	slope_limit = util.GetParamNumber("-slope_limit", 2e-02, "regularization factor in slip and diff velocity"),
	granular_model= util.GetParamNumber("-granular_model", 3, "Opt: 0 Const, 1 Linear, 2 Einstein, 3 Rheology(I) + Einstein, 4 Relax"),
	density_model  = util.GetParam("-density_model", "linear", "constant, linear"),
	drag_mod = util.GetParamNumber("-drag_model", 2, "Opt: 0 StokesLaw, 1 formula, 2 Schiller-Naumann, 3 Turton and Levenspiel"),
	riemman = util.GetParamNumber("-riemman", 1, "Opt: 0 Godunov, 1 Rusanov, 2 Roe"),
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
	roughness_length = 1e-04,
	
}

params.startTime  = 0.0
params.endTime    = params.DT * params.numTimeSteps
params.DTmax = params.DT
params.DTLimex = params.DT

c_init = params.c_init
params.interface_value  = params.alpha_min/params.packing_factor


------------------------------------------------------------------------------------------
-- Geometry Parameters  - GridName -  Domain Subsets
------------------------------------------------------------------------------------------

local geometry = myProblem:GetGeometry(params)
local requestedGrid = util.GetParam("-geom", "")
if requestedGrid ~= "" then geometry.gridName = requestedGrid end
params.gridName = geometry.gridName


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

--InitUG (params.dim, AlgebraType("CPU", params.dim+2))
InitUG (params.dim, AlgebraType("CPU", 1))


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

	approxSpace,u = myProblem:ApproximationSpace(geometry.allSubsets)
	
------------------------------------------------------------------------------------------
-- Lua Functions
------------------------------------------------------------------------------------------

myProblem:RegisterCallbacks()

-------------------------------------------------------------------------- Parameters List
------------------------------------------------------------------------------------------
-- Parameters List
------------------------------------------------------------------------------------------


InterfaceValues = myProblem:InterfaceParameters()


------------------------------------------------------------------------------------------
-- Secondary Variables (Closures)
------------------------------------------------------------------------------------------


myProblem:Clousures(approxSpace,u,geometry.turbulenceZeroSubsets)


------------------------------------------------------------------------------------------
-- Compose the discretization
------------------------------------------------------------------------------------------


NavierStokesDisc = myProblem:Discretization(geometry.innerSubsets)


------------------------------------------------------------------------------------------
-- Boundary Conditions
------------------------------------------------------------------------------------------


local boundaries = myProblem:CreateBoundaryConditions(NavierStokesDisc, geometry)


---------------------------------------------------------------------------------------
-- Parameters Inputs
---------------------------------------------------------------------------------------

myProblem:ConnectClosures(NavierStokesDisc)

---------------------------------------------------------------------------------------
-- Global Discretization
---------------------------------------------------------------------------------------

local domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(NavierStokesDisc)
for _, boundary in ipairs(boundaries) do domainDisc:add(boundary) end

---------------------------------------------------------------------------------------
-- Time Discretization
---------------------------------------------------------------------------------------

print("Setting Time Discretization")

local timeDisc = myProblem:TimeDiscretization(domainDisc)


------------------------------------------------------------------------------------------
-- Interpolate initial values
------------------------------------------------------------------------------------------
print("Initializing Values")

local time, step, time_work_total = myProblem:InitializeSolution(u, folder_vtk)

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
if params.boolAveNormal then
	myProblem.Normal:update()
end

------------------------------------------------------------------------------------------
-- Steady State Solution
------------------------------------------------------------------------------------------
print("Calculating SteadyState")
local Newton_Steps = 0
local Newton_Steps_fail = 0
local time_work_steady=0.0
local linsolver_calls = 0
local linsolver_steps = 0

if params.doSteadyState then
	-- Steady state solution.
	
	NewtonSolverSteady:add_inner_step_update(myProblem.gamma)
	NewtonSolverSteady:add_step_update(myProblem.RelVel)
	if params.boolAveNormal then
		NewtonSolverSteady:add_step_update(myProblem.Normal)
	end
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
	Newton_Steps = 1
	Newton_Steps_fail = Newton_Steps - boolSolution
end

if(params.boolFixVel) then
	fixer = DirichletBoundary()
	domainDisc:add(fixer)
	fixer:invert_subset_selection()
	fixer:add("u", "")
	fixer:add("v", "")
	fixer:add("p", "")
	if params.dim == 3 then
		fixer:add("w", "")
	end
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
if (boolAveNormal) then
	NLSolver:add_step_update(myProblem.Normal)
end

if params.timeMethod == "limex" then
	NLSolver:add_step_update(myProblem.gamma)
	if params.boolSlipDiff then
		NLSolver:add_step_update(myProblem.SlipDiff)
	elseif params.boolSlipVel then
		NLSolver:add_step_update(myProblem.SlipVel)
	end
else
	NLSolver:add_inner_step_update(myProblem.gamma)
	if (boolAveNormal) then
		NLSolver:add_inner_step_update(myProblem.Normal)
	end
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
	if (step % params.outputFactor == 0 ) then
		local vtkStep = math.floor(step / params.outputFactor)
		print("Writing initial values")
		out:print_subsets(vtk_file_name, u,geometry.allSubsets,vtkStep,time, true)
		print ("Output to file " .. vtk_file_name .. ".vtu  in time t =" .. time)
		print ("    -   -   -   -   -   -   -   -   -   -   -   -   -   -   ")
		print ("                                                            ")
		print ("    -   -   -   -   -   -   -   -   -   -   -   -   -   -   ")
	end
	
	myProblem:SaveCheckPoint(u,folder_vtk)
	if  (params.writeIntegral and step==0) then
	
		local Value_inner1, Value_inner2 = myProblem:ComputeIntegrals(u, geometry.innerSubsets)
		
		if(rank == 0) then
			myProblem:WriteValues( folder_vtk, step, time, Value_inner1, Value_inner2, time_work_steady, time_work_total, Newton_Steps, Newton_Steps_fail, linsolver_calls, linsolver_steps,false)
		end
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
tBefore = os.clock()

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
			Newton_Steps, Newton_Steps_fail, linsolver_calls_step, linsolver_steps_step , boolSolution = myProblem:SolveNonlinearProblem( u, NLSolver, op, solTimeSeries, params.DT, step,StartTime,EndTime)
		end
		time = EndTime
		tAfter_step = os.clock()
		
		
		if boolSolution == 1 then
		
			if (step % params.outputFactor == 0 ) then
				local vtkStep = math.floor(step / params.outputFactor)
				out:print_subsets(vtk_file_name, u,geometry.allSubsets,vtkStep,time)
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
					

			
			
			
			if (params.writeIntegral) then
			
				local Value_inner1, Value_inner2 = myProblem:ComputeIntegrals(u, geometry.innerSubsets)
				time_work_total = time_work_total + tAfter_step - tBefore_step
				if(rank==0) then
					myProblem:WriteValues( folder_vtk, step, time, Value_inner1, Value_inner2, tAfter_step - tBefore_step, time_work_total, Newton_Steps, Newton_Steps_fail, linsolver_calls_step, linsolver_steps_step,false)
				end
			end
						
		else
			print("++++++ TIMESTEP " .. step .. "  FAILED ++++++")
			local vtkStep = math.floor(step / params.outputFactor)
			out:print_subsets(vtk_file_name, u,geometry.allSubsets,vtkStep,time)
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
	print("Total Computation took " .. time_work_total .. " seconds.")
	print("-------------------------------------------------------------------------------")
	print("")
	print("")
	print ("Output to file " .. vtk_file_name .. ".vtu")
	print("done.")


	if (params.writeIntegral) then
	
		local Value_inner1, Value_inner2 = myProblem:ComputeIntegrals(u, geometry.innerSubsets)
		if(rank == 0) then
			myProblem:WriteValues( folder_vtk, params.numTimeSteps, time, Value_inner1, Value_inner2, time_work_total, time_work_total, total_Newton_Steps, total_Newton_Steps_fail, total_linsolver_calls_step, total_linsolver_steps_step,true)
		end
	end
end

SynchronizeProcesses()
if (params.NewtonDebug and rank == 0 and SpaceSize > 1) then

	csvfile.MergeDebugPVD(myProblem.debug_dir, params.file_name)
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

