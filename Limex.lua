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
local numPreRefs = util.GetParamNumber("--numPreRefs", 0, "number of refinements before parallel distribution")
local numRefs    = util.GetParamNumber("--num-refs",    2, "number of refinements")
local startTime  = util.GetParamNumber("--start", 0.0, "start time")
local endTime    = util.GetParamNumber("--end", 1.0, "end time")
local numTimeSteps    = util.GetParamNumber("--n_steps", 5, "time steps")

local doVTK      = util.GetParamNumber("--with-vtk",false,"Save Values")

params.tol     = util.GetParamNumber("--limex-tol", 1e-2, "time step size")
params.nstages = util.GetParamNumber("--limex-nstages", 2, "limex stages (2 default)")
params.limex_partial_mask = util.GetParamNumber("--limex-partial", 0, "limex partial (0 or 3)")
params.limex_debug_level = util.GetParamNumber("--limex-debug-level", 5, "limex debug level (integer)")


dt=(endTime-startTime)/numTimeSteps
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
--if dim == 2 then gridName = util.GetParam("-grid", "unit_square_01_quads_2x2.ugx")
if dim == 2 then gridName = util.GetParam("-grid", "Dune2D_tri_double.ugx")
else print("Dimension "..dim.." not supported."); exit(); end

print(" Selected Parameter:")
print("    numRefs      = " .. numRefs)
print("    numPreRefs   = " .. numPreRefs)
print("    startTime 	  = " .. startTime)
print("    endTime 		  = " .. endTime)
print("    dt 			    = " .. dt)
print("    eps          = " .. eps)
print("    grid         = " .. gridName)


print("    limex_nstages      = " .. params.nstages)
print("    limex_partial_mask = " .. params.limex_partial_mask)
print("    limex_tol          = " .. params.tol)

-- choose algebra
InitUG(dim, AlgebraType("CPU", 1));

-- Create, Load, Refine and Distribute Domain
local mandatorySubsets = {"Inner","Inner2", "Left" , "Top", "Bottom", "Right"}
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
approxSpace:add_fct("c", "Lagrange", 1)

-- lets order indices using Cuthill-McKee
--OrderCuthillMcKee(approxSpace, true);

--------------------------------------------------------------------------------
--  Setup FV Convection-Diffusion Element Discretization
--------------------------------------------------------------------------------

print (">> Setting up Assembling")


-- The coordinates (cx, cy) specify the rotation center of the cone
local cx = 20.5
local cy = 7.5

-- The coordinates (ax, ay) specify the position of the highest point of the
-- cone at start time t=0.0
local ax = 0.25
local ay = 0.0

-- The parameter nu specifies the rotation velocity
local nu = 1.0 -- 100

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
	return	nu*(y - cx), nu*(cy - x)
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
     vtk:print("ConvDiffSol", u, step, time)
     --local ref = u:clone()
     --Interpolate("exactSolution", ref, "c", time)
     --vtk:print("ConvDiffRef", ref, step, time)
  end
end



-- grid function debug writer
local dbgWriter = GridFunctionDebugWriter(approxSpace)
GetLogAssistant():set_debug_level("LIB_LIMEX", 4)
-- descriptor for linear solver
local solverDesc = {
    name = "bicgstab", -- "linear"
    precond = {
      type = "gmg",
      approxSpace = approxSpace,
      smoother = "ilu",
      rap=true,

      preSmooth=2,
      postSmooth=2,
      baseLevel = 2,
      baseSolver = "lu",
    },
    convCheck = {
      type ="standard",
      maxSteps = 100,
      minDef = 1e-9,
      reduction = 1e-12 },
      
      
}


print (">> Setting up Assembling (2)")
--------------------------------------------------------------------------------
--  Standard schemes
--------------------------------------------------------------------------------
local elemDisc ={}
local dirichletBND = {}
local domainDisc = {}


local upwind = FullUpwind()
--local upwind = NoUpwind()
-- setup for discretizations

elemDisc = ConvectionDiffusion("c", "Inner,Inner2", "fv1")
elemDisc:set_upwind(upwind)
elemDisc:set_diffusion(eps)
elemDisc:set_velocity("Velocity")


dirichletBND = DirichletBoundary()
dirichletBND:add(0.0, "c", "Bottom,Left")
dirichletBND:add(0.5, "c", "Top,Right")

domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(elemDisc)
domainDisc:add(dirichletBND)



local limexLSolver = nil
local limexNLSolver = nil

local limexConvCheck=ConvCheck(1, 5e-8, 1e-10, true)
limexConvCheck:set_supress_unsuccessful(true)


limexLSolver = util.solver.CreateSolver(solverDesc)
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

local tSteps = (endTime-startTime)/numTimeSteps
local dtlimex = tSteps

local gridSize = 0.5*math.pow(0.5, numRefs)
if (doVelocity) then
  local tCFL = gridSize/nu
  dtlimex = math.min(dtlimex, tCFL)
end
--  Euclidean (algebraic) norm
--local estimator = Norm2Estimator()
--tol = 0.37/(gridSize)*tol


--print (estimator)
local limexEstimator = CompositeGridFunctionEstimator()
limexEstimator:add(L2ComponentSpace("c", 2))

-- descriptor for integrator
local limexDesc = {

  nstages = params.nstages,
  steps = {1,2,3,4,5,6},
  domainDisc=domainDisc,
  nonlinSolver = limexNLSolver,
  
  tol = params.tol,
  dt = dtlimex,
  dtmin = 1e-9,
  
  rhoSafetyOPT = 0.25,
  
}


-- setup for time integrator
local limex = util.limex.CreateIntegrator(limexDesc)

limex:set_dt_min(1e-9)
limex:add_error_estimator(limexEstimator)
limex:set_increase_factor(2.0)

if (vtk) then
	
	--limex:attach_observer(vtkObserver)
	limex:attach_observer(luaObserver)
	print("-------------------------------------------------------------")
	print("-------------------------------------------------------------")
	print("-------------------------------------------------------------")
	print("-------------------------------------------------------------")
         
end


--limex:attach_observer(refObserver)



limex:set_stepsize_greedy_order_factor(1.0)
limex:select_cost_strategy(LimexNonlinearCost())
-- limex:disable_matrix_cache()  -- recompute ()
limex:enable_matrix_cache() -- keep matrix

print ("dtLimex   = "..dtlimex)
print ("hGrid     = "..gridSize)
print ("tolLimex  = "..params.tol)

dbgWriter:set_vtk_output(false)
dbgWriter:set_conn_viewer_output(true)
-- limexNLSolver:set_debug(dbgWriter)

limex:set_stepsize_greedy_order_factor(1.0)
limex:select_cost_strategy(LimexNonlinearCost())
-- limex:disable_matrix_cache()  -- recompute ()
limex:enable_matrix_cache() -- keep matrix

-- set initial value
print(">> Interpolating start values")
local u = GridFunction(approxSpace)
u:set(0.0)
Interpolate("exactSolution", u, "c", startTime)
-- solve problem

print(">> Peclet number:"..nu*1.0/eps)
print(">> Grid Peclet number:"..nu*gridSize/eps)
print(">> Solve problem")

out = VTKOutput()
out:clear_selection()
out:select_all(false)
out:select_nodal("c", "c")
out:print("Sol", u,0,startTime)


local cstart=os.clock()
for step = 1, numTimeSteps do
	StepTime=startTime+(step)*dt
	limex:apply(u, StepTime, u, StepTime-dt)
	out:print("Sol", u,step,StepTime)
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
