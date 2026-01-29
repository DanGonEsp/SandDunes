--------------------------------------------------------------------------------
--[[!
-- \file sample/conv-diff-limex.lua
-- \ingroup app_convdiff
-- \{
-- \author Arne Naegel
-- \brief Illustrates limex usage based on conv_diff.lua by Andreas Vogel
-- \}
-- ]]
--------------------------------------------------------------------------------

ug_load_script("ug_util.lua")
ug_load_script("plugins/Limex/lua/limex_util.lua")
ug_load_script("util/profiler_util.lua")
ug_load_script("util/load_balancing_util.lua")

local params = {}
local dim = util.GetParamNumber("-dim", 2, "dimension")
local file_name = util.GetParam("-file_name", "Limex")
local numPreRefs = util.GetParamNumber("--numPreRefs", 0, "number of refinements before parallel distribution")
local numRefs    = util.GetParamNumber("--num-refs",    1, "number of refinements")
local startTime  = util.GetParamNumber("--start", 0.0, "start time")
local endTime    = util.GetParamNumber("--end", 10000000.0, "end time")
local numTimeSteps    = util.GetParamNumber("--n_steps", 5, "time steps")

local doVTK      = util.GetParamNumber("--with-vtk",true,"Save Values")

StatBool = util.GetParamBool("-StatBool", true)
max_newton_steps_steady_state=util.GetParamNumber("-numNewtonSteps", 100)
max_newton_steps_transcient=util.GetParamNumber("-max_newton_steps_transcient", 50)
max_linear_steps=util.GetParamNumber("-max_linear_steps", 200)
AbsDefect = util.GetParamNumber("-AbsDefect", 1e-05)
RedDefect = util.GetParamNumber("-RedDefect", 1e-03)
damping_mg = util.GetParamNumber("-damping_mg", 0.5)
value_beta = util.GetParamNumber("-value_beta", -0.01)
lambdamaxSteps = util.GetParamNumber("-lambdamaxSteps", 7)
lambdaStart  = util.GetParamNumber("-lambdaStart", 1.0)

params.tol     = util.GetParamNumber("--limex-tol", 1e-2, "time step size")
params.nstages = util.GetParamNumber("--limex-nstages", 2, "limex stages (2 default)")
params.limex_partial_mask = util.GetParamNumber("--limex-partial", 0, "limex partial (0 or 3)")
params.limex_debug_level = util.GetParamNumber("--limex-debug-level", 5, "limex debug level (integer)")


-- This scales the amount of diffusion of the problem
local eps       = util.GetParamNumber("--eps", 1e-1, "strength of diffusion")
local doVelocity = true

-- if withRedist == true, you may alter the redistribution behavior over the
-- parameters declared in util/load_balancing_util.lua.
-- Note that numPreRefs is ignored in this case.
local withRedist = util.HasParamOption("--withRedist", false, "Enables a more sophisticated distribution approach.")

GetLogAssistant():set_debug_level("LIB_LIMEX", params.limex_debug_level)


util.CheckAndPrintHelp("Time-dependent problem setup example\n(by Andreas Vogel)");

local gridName = nil
if dim == 2 then gridName = util.GetParam("-grid", "Dune2D_tri_double.ugx")
else print("Dimension "..dim.." not supported."); exit(); end

------------------------------------------------------------------------------------------
-- Folder and files
------------------------------------------------------------------------------------------

local folder = file_name
local vtk_file_name = file_name
if not DirectoryExists (folder) then
    CreateDirectory (folder)
end
vtk_file_name = folder .. "/" .. vtk_file_name -- VTK output file name base

print(" Selected Parameter:")
print("    numRefs      = " .. numRefs)
print("    numPreRefs   = " .. numPreRefs)
print("    startTime 	  = " .. startTime)
print("    endTime 		  = " .. endTime)
print("    numTimeSteps   = " .. numTimeSteps)
print("    eps          = " .. eps)
print("    grid         = " .. gridName)


print("    limex_nstages      = " .. params.nstages)
print("    limex_partial_mask = " .. params.limex_partial_mask)
print("    limex_tol          = " .. params.tol)

-- choose algebra
InitUG(dim, AlgebraType("CPU", 1));

if dim == 3 then
	fct_cmp_tbl = {"u", "v", "w", "p", "c"}
	vel_cmp_tbl = {"u", "v", "w"}
else if dim == 2 then
		fct_cmp_tbl = {"u", "v", "p", "c"}
		vel_cmp_tbl = {"u", "v"}
	else print("Dimension "..dim.." not supported."); exit(); end
end
	

-- Create, Load, Refine and Distribute Domain
local mandatorySubsets = {"Inner","Inner2", "Left" , "Top", "Bottom", "Right"}
local allSubsets = "Inner,Inner2,Left,Top,Bottom,Right"
local Inner_total={"Inner","Inner2"}
local dom
if withRedist == true then
	dom = util.CreateDomain(gridName, 0, mandatorySubsets)
	balancer.RefineAndRebalanceDomain(dom, numRefs)
else
	dom = util.CreateAndDistributeDomain(gridName, numRefs, numPreRefs, mandatorySubsets)
end

print("\ndomain-info:")
print(dom:domain_info():to_string())

-- create Approximation Space
print(">> Create ApproximationSpace")
local approxSpace = ApproximationSpace(dom)

approxSpace:add_fct("u", "Lagrange",1,allSubsets)
approxSpace:add_fct("v", "Lagrange",1,allSubsets)
if dim == 3 then
	approxSpace:add_fct("w", "Lagrange",1,allSubsets)
end
approxSpace:add_fct("p", "Lagrange",1,allSubsets)
approxSpace:add_fct("c", "Lagrange",1,allSubsets)

-- lets order indices using Cuthill-McKee
--OrderCuthillMcKee(approxSpace, true);

--------------------------------------------------------------------------------
--  Setup FV Convection-Diffusion Element Discretization
--------------------------------------------------------------------------------

print (">> Setting up Assembling")


-- The coordinates (cx, cy) specify the rotation center of the cone
local cx = 0.0
local cy = 7.5

-- The coordinates (ax, ay) specify the position of the highest point of the
-- cone at start time t=0.0
local ax = 0.25
local ay = 0.0

-- The parameter nu specifies the rotation velocity
local nu = 10.0 -- 100

-- The parameter delta is a scaling factor influencing the steepness of the cone
delta = 1e+2  --1e-1
 
-- This is the exact solution for our problem
function exactSolution(x, y, t)
	local xRot = math.cos(nu*t) * (x-cx) - math.sin(nu*t) * (y-cy) 
	local yRot = math.sin(nu*t) * (x-cx) + math.cos(nu*t) * (y-cy) 
	
	local expo = -((xRot - ax)*(xRot - ax) + (yRot - ay)*(yRot - ay)) / (delta + 4*eps*t)
	local scale = delta/(delta+4*eps*t)

	return scale * math.exp(expo)
end
	
-- The velocity field
function Velocity(x, y, t)
	--return	nu*(y - cx), nu*(cy - x)
	return	nu, 0.0
end
	
-- The dirichlet condition
function DirichletValue(x, y, t, si)
	return true, exactSolution(x, y, t)
end


local vtk = nil 
local vtkObserver = nil
 
if (doVTK) then
   print("doVTK = true")
   vtk = VTKOutput();
   vtkObserver=VTKOutputObserver("MyFile.vtk", vtk)
end

-- post-processing (after each step)
function postProcess(u, step, time, currdt)
  local l2error = L2Error("exactSolution", u, "c", time, 4)
  local l2norm=L2Norm(u,"c", time)
  print("L2Error at t=\t"..time.."\t=\t"..l2error..", (rel=".. (l2error/l2norm).. ")\tL2Norm=\t"..l2norm)

  if (doVTK) then
     --vtk:print("ConvDiffSol", u, step, time)
     --local ref = u:clone()
     --Interpolate("exactSolution", ref, "c", time)
     --vtk:print("ConvDiffRef", ref, step, time)
  end
end

---------------------------------------------------------------------- Interface Properties

InterfaceValues = Interface()
InterfaceValues:set_bool_initialized(true)

print (">> Setting up Assembling (2)")
--------------------------------------------------------------------------------
--  Standard schemes
--------------------------------------------------------------------------------
local elemDisc ={}
local dirichletBND = {}
local domainDisc = {}

-- setup for discretizations

elemDisc =  NavierStokesFV1M (fct_cmp_tbl, Inner_total)
elemDisc:set_exact_jacobian (false)
elemDisc:set_stokes (false)
elemDisc:set_laplace ( false)
elemDisc:set_upwind ("full")
elemDisc:set_upwind_vol("full")
elemDisc:set_peclet_blend (false)
elemDisc:set_stabilization ("fields_2", "raw")

elemDisc:set_density(1000)
elemDisc:set_kinematic_viscosity (1e-06)

elemDisc:set_phase_parameters(InterfaceValues)




InletDisc = NavierStokesInflowFV1M (elemDisc)
InletDisc:add ("Velocity","Left")

dirichletBND = DirichletBoundary()
dirichletBND:add(0.0, "c", "Left")
dirichletBND:add(1.0, "c", "Bottom,Top")


-- boundary condition at the outlet
OutletDisc = NavierStokesNoNormalStressOutflowFV1M (elemDisc)
OutletDisc:add ("Right")
OutletDisc:set_phase_parameters(InterfaceValues)

-- boundary condition at the impermeable walls
WallDisc = NavierStokesWall (elemDisc)
WallDisc:add ("Bottom,Top")


domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(elemDisc)
domainDisc:add (InletDisc)
domainDisc:add (OutletDisc)
domainDisc:add (WallDisc)
domainDisc:add(dirichletBND)




------------------------------------------------------------------------------------------
-- Set up the solver
------------------------------------------------------------------------------------------

-- grid function debug writer
local dbgWriter = GridFunctionDebugWriter(approxSpace)
GetLogAssistant():set_debug_level("LIB_LIMEX", 4)

-- descriptor for linear solver
local LinearSolverDesc = {
	name = "bicgstab",
	precond =
	{
		type = "gmg",
		approxSpace = approxSpace,
		rap=true,
		rim = false,
		cycle = "V",
		smoother =
		{
			type = "ilu",
			beta = value_beta,
			damping 	= damping_mg,
			sort	= false,
			--sortEps 	= 1.e-50,
			inversionEps 	= 1.e-16,
			consistentInterfaces   = false,     --consistentInterfaces and overlap shouldnot be activated at the same time
			overlap 		= true,             --consistentInterfaces and overlap shouldnot be activated at the same time
			--ordering 		= nil
		},

		preSmooth=2,
		postSmooth=2,
		baseLevel = numPreRefs,
		baseSolver = "lu",
    },
	convCheck =
	{
		type		= "standard",
		iterations	= max_linear_steps,
		absolute	= 1e-12,
		reduction	= 1e-3,
		verbose		= true
	}
      
      
}
local NewtonSolverDesc =
{
	type = "newton",
	debug = false, -- for the debug output from the Newton's method
	linSolver = LinearSolverDesc,
	lineSearch =
	{
		type			= "standard",
		maxSteps		=lambdamaxSteps,
		lambdaStart		= lambdaStart,
		lambdaReduce	= 0.7,
		acceptBest 		= true,
		checkAll		= false,
		suffDesc		= 0.3,
		maxDefect	= 2e20
		
	},
	convCheck =
	{
		type		= "standard",
		iterations	= max_newton_steps_steady_state,
		absolute	= AbsDefect,
		reduction	= RedDefect,
		verbose		= true
	}
}

local NewtonSolverSteady = util.solver.CreateSolver(NewtonSolverDesc)
local limexLSolver = util.solver.CreateSolver(LinearSolverDesc)


local limexConvCheck=ConvCheck(1, 5e-5, 1e-8, true)
limexConvCheck:set_supress_unsuccessful(true)


local limexNLSolver = nil
limexNLSolver = NewtonSolver()
limexNLSolver:set_linear_solver(limexLSolver)
limexNLSolver:set_convergence_check(limexConvCheck)
print(limexNLSolver)



-- local refObserver = PlotRefOutputObserver("DirichletValue", vtk) -- now obsolete
local luaObserver = LuaCallbackObserver()

-- work-around (waiting for implementation of SmartPtr forward to lua...)
function luaPostProcess(step, time, currdt)
  print("LUAPostProcess: "..step..","..time..","..currdt)
  postProcess(luaObserver:get_current_solution(), step, time, currdt)
  return 0;
end
luaObserver:set_callback("luaPostProcess")

local dtlimex = (endTime-startTime)/numTimeSteps

local gridSize = 1.0
--  Euclidean (algebraic) norm
--local estimator = Norm2Estimator()
--tol = 0.37/(gridSize)*tol


--print (estimator)
local limexEstimator = CompositeGridFunctionEstimator()
--limexEstimator:add(H1SemiComponentSpace("u", 2 ))
--limexEstimator:add(H1SemiComponentSpace("v", 2 ))

--limexEstimator:add(L2ComponentSpace("u", 2))
--limexEstimator:add(L2ComponentSpace("v", 2))

--limexEstimator:add(H1SemiComponentSpace("p", 2))
limexEstimator:add(L2ComponentSpace("c", 2))

-- descriptor for integrator
local limexDesc = {

  nstages = params.nstages,
  steps = {1,2,3,4,5,6},
  nthreads = 1,
  domainDisc=domainDisc,
  nonlinSolver = limexNLSolver,
  -- makeConsistent = true,
  
  tol = params.tol,
  dt = dtlimex,
  dtmax = dtlimex,
  dtmin = 1e-5,
  rhoSafetyOPT = 0.25,
  dtred = 0.5,
  dtIncr = 2.0,
  matrixCache = true,
  conservative = false
  
}


-- setup for time integrator
local limex = util.limex.CreateIntegrator(limexDesc)

limex:set_dt_min(limexDesc.dtmin)
limex:set_dt_max(limexDesc.dtmax)
limex:set_reduction_factor(limexDesc.dtred)
limex:set_increase_factor(limexDesc.dtIncr)
limex:add_error_estimator(limexEstimator)
limex:set_stepsize_greedy_order_factor(1.0)
limex:select_cost_strategy(LimexNonlinearCost())


if (false) then
	
	--limex:attach_observer(vtkObserver)
	limex:attach_observer(luaObserver)
end


--limex:attach_observer(refObserver)




print ("dtLimex   = "..dtlimex)
print ("hGrid     = "..gridSize)
print ("tolLimex  = "..params.tol)

dbgWriter:set_vtk_output(false)
dbgWriter:set_conn_viewer_output(true)
-- limexNLSolver:set_debug(dbgWriter)

-- set initial value
print(">> Interpolating start values")
local u = GridFunction(approxSpace)
u:set(0.0)
--Interpolate("exactSolution", u, "c", startTime)
-- solve problem

print(">> Peclet number:"..nu*1.0/eps)
print(">> Grid Peclet number:"..nu*gridSize/eps)
print(">> Solve problem")
	
if StatBool then
	-- Fix the mass fraction and solve the linear problem for the momentum
	fixer = DirichletBoundary()
	domainDisc:add(fixer)
	fixer:invert_subset_selection()
	fixer:add("c", "")

	NewtonSolverSteady:init(AssembledOperator(domainDisc))


	
	NewtonSolverSteady:prepare(u)

	-- apply the solver for the stationary pressure problem
     print("++++++ STEADY STATE CALCULATION BEGIN ++++++")
	tBefore_s= os.clock()
	if not NewtonSolverSteady:apply(u) then
		print("===> THE PREPARATION PHASE FAILED! <===")
		exit()
	end
	tAfter_s = os.clock()
    num_newton_steps = NewtonSolverSteady:num_newton_steps()
    linsolver_calls = NewtonSolverSteady:total_linsolver_calls()
    linsolver_steps = NewtonSolverSteady:total_linsolver_steps()
    average_linear_steps = NewtonSolverSteady:total_average_linear_steps()
    average_non_linear_rates = NewtonSolverSteady:total_average_non_linear_rates()
    
    print("num_newton_steps = " .. num_newton_steps .. ".")
    print("linsolver_calls = " .. linsolver_calls .. ".")
    print("linsolver_steps = " .. linsolver_steps .. ".")
    print("average_linear_steps = " .. average_linear_steps .. ".")
    print("average_non_linear_rates = " .. average_non_linear_rates .. ".")
    
    NewtonSolverSteady:clear_average_convergence();
    
    time_work_step = tAfter_s-tBefore_s
	print("Computation for steady state took " .. time_work_step .. " seconds.")
	domainDisc:remove (fixer)
	
	print("++++++++++++++++++++++++ INITIAL CONDITIONS  (STEADY STATE DONE) ++++++++++++++++++++++++")
	

end




if (doVTK) then
	out = VTKOutput()
	out:clear_selection()
	out:select_all(false)
	out:select_nodal("u", "u")
	out:select_nodal("v", "v")
	out:select_nodal("p", "p")
	out:select_nodal("c", "c")
	out:print(vtk_file_name, u,0,startTime)
end


local cstart=os.clock()
for step = 1, numTimeSteps do
	print("---------------------------------------------								++++++ TIMESTEP " .. step .. " BEGIN ++++++")

	
	StepTime=startTime+(step)*dtlimex
	limex:apply(u, StepTime, u, StepTime-dtlimex)
	if (doVTK) then
		out:print(vtk_file_name, u,step,StepTime)
	end
	
	print("---------------------------------------------								++++++ TIMESTEP " .. step .. "  END ++++++	TIME =  "..StepTime)

end
local cend=os.clock()


print ("CDELTA=\t"..cend - cstart)



print("Writing profile data")
WriteProfileData("profile_data.pdxml")
util.PrintProfile_TotalTime("main ")

FreeUserData()

-- end group app_convdiff
--[[!
\}
]]--
