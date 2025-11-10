--------------------------------------------------------------------------------
--
--   Lua - Script to compute the cylinder problem
--
--	This script sets up a problem for the Navier-Stokes discretization
--	and solves the cylinder problem
--
--   Author: Arne Naegel
--
--------------------------------------------------------------------------------

local myPath = ug_get_current_path()
package.path = package.path..";".. myPath.."/?.lua"

ug_load_script("ug_util.lua")
ug_load_script ("util/load_balancing_util.lua")
-- ug_load_script("util/domain_disc_util.lua")
-- ug_load_script("util/conv_rates_static.lua")
---ug_load_script("plugins/Limex/limex_util.lua")
ug_load_script("navier_stokes_util.lua")
ug_load_script ("ug_util.lua")
ug_load_script ("util/load_balancing_util.lua")
RequiredPlugins({"Limex", "NavierStokes"})



-- Command line arguments
local ARGS = {
  dim     = util.GetParamNumber("-dim", 2, "world dimension"),
  doSteadyState = util.GetParam("--steady-state",false, "Compute steady state solution"),
  doInitialSteadyStateValues = util.GetParam("--steadyInitialValues",true, "Compute steady state solution"),
  bStokes   = util.HasParamOption("-stokes", "If defined, only Stokes Eq. computed"),
  bExactJac   = util.HasParamOption("-exactjac", "If defined, exact jacobian used"),
  bLaplace  = not util.HasParamOption("-nolaplace", "If defined, only laplace term used"),
  bPecletBlend= util.HasParamOption("-pecletblend", "If defined, Peclet Blend used"),
  
  upwind      = util.GetParam("-upwind", "full", "Upwind type"),
  
  stab        = util.GetParam("-stab", "fields_2", "Stabilization type"),
  diffLength  = util.GetParam("-difflength", "cor", "Diffusion length type"),
  
  stabGrad       = util.GetParamNumber("--stabGrad", 0.1, "Stabilization parameter."),
  stabStreamline = util.GetParamNumber("--stabStreamline", 0.0, "Stabilization parameter."),
  stabDiv        = util.GetParamNumber("--stabDiv", 0.0, "Stabilization parameter."),
 

  limexNStages = util.GetParamNumber("--limex-num-stages", 4, "number of LIMEX stages: 0:steady state, 1:std. integration, <= 2:LIMEX"),
  limexDebugLevel = util.GetParamNumber("--limex-debug-level", 2, "debug level"),
  limexTOL = util.GetParamNumber("--limex-tol", 1e-2, "debug level"),
  
  solverID =  util.GetParam("--solver-id", "gmg", "superlu, gmg"),
  
  numRefs   = util.GetParamNumber("--numRefs", 2, "number of grid refinements"),
  numPreRefs  = util.GetParamNumber("-numPreRefs", 0, "number of prerefinements (parallel)")
  
 
  
}

ARGS.discType  = util.GetParam("-type", "fv1", "Disc type 		fv1	fv	fe	fecr	 fvcr")
ARGS.vorder, ARGS.porder = 1,1





-- Shortcut.
local dim =ARGS.dim

if dim == 2 then 
	ARGS.gridName = util.GetParam("-grid", "grids/cylinderp.ugx")
elseif dim == 3 then
	ARGS.gridName = util.GetParam("-grid", "grids/cylinder3d.ugx")
else print("Selected Dimension not supported. Exiting."); exit(); end


------------------------------------------------------------------------------------------
-- Get command line parameters
------------------------------------------------------------------------------------------

-- Physical parameters
viscosity 	= util.GetParamNumber("-visc", 1e-3, "kinematic viscosity")
--inflow		= util.GetParamNumber("-inflow", 10, "max. inflow velocity")
bStokes 	= util.HasParamOption("-Stokes", "If defined, only Stokes Eq. computed")

-- Numerical parameters of the discretization
numRefs 	= util.GetParamNumber("-numRefs", 2, "number of grid refinements")
numPreRefs 	= util.GetParamNumber("-numPreRefs", 0, "number of prerefinements (parallel)")
bNoLaplace 	= util.HasParamOption("-noLaplace", "If defined, only laplace term used")
bExactJac 	= util.HasParamOption("-exactJac", "If defined, exact jacobian used")
bPecletBlend= util.HasParamOption("-PecletBlend", "If defined, Peclet Blend used")
upwind      = util.GetParam("-upwind", "full", "Upwind type (no, full, weighted, lps, pos, reg)")
bPac        = util.HasParamOption("-pac", "If defined, pac upwind used")
stab        = util.GetParam("-stab", "fields_2", "Stabilization type (fields or flow)")
diffLength  = util.GetParam("-difflength", "cor", "Diffusion length type (raw, fivepoint or cor)")

-- Parameters of the solver
ilu_beta	= util.GetParamNumber("-iluBeta", -0.8, "choose a negative value depending on the convection rate")

------------------------------------------------------------------------------------------
-- Geometry data
------------------------------------------------------------------------------------------

dim 		= 2
gridName = "cylinderp.ugx"
inflow_y_0	= 0 -- y coord. of the lower corner of the inflow
inflow_y_1	= 0.41 -- y coord. of the upper corner of the inflow
inflow = 10


------------------------------------------------------------------------------------------
-- Print the parameters
------------------------------------------------------------------------------------------

-- Lets write some info about the choosen parameter
print(" Selected Parameter:")
print("    dim              = " .. dim)

print("    numTotalRefs     = " .. ARGS.numRefs)
print("    numPreRefs       = " .. ARGS.numPreRefs)
print("    grid             = " .. ARGS.gridName)

print("    porder           = " .. ARGS.porder)
print("    vorder           = " .. ARGS.vorder)
print("    type             = " .. ARGS.discType)
print("    only stokes      = " .. tostring(ARGS.bStokes))
print("    only laplace     = " .. tostring(not(ARGS.bLaplace)))
print("    exact jacobian   = " .. tostring(ARGS.bExactJac))
print("    peclet blend     = " .. tostring(ARGS.bPecletBlend))
print("    upwind           = " .. ARGS.upwind)
print("    stab             = " .. ARGS.stab)
print("    stabGrad         = " .. ARGS.stabGrad)
print("    stabDiv          = " .. ARGS.stabDiv)
print("    stabStreamline   = " .. ARGS.stabStreamline)
print("    diffLength       = " .. ARGS.diffLength)



print("    SteadyState      = " .. tostring(ARGS.doSteadyState))
print("    LIMEX.TOL        = " .. ARGS.limexTOL)
print("    LIMEX.NumStages  = " .. ARGS.limexNStages)

--------------------------------------------------------------------------------
-- Debug output.
--------------------------------------------------------------------------------
--local logAssistant = GetLogAssistant()
--logAssistant:set_debug_level("LIB_LIMEX", ARGS.limexDebugLevel)
--logAssistant:set_debug_levels(10)
--logAssistant:set_debug_level("LIB_DISC_MULTIGRID", 10)

--------------------------------------------------------------------------------
-- InitUG
--------------------------------------------------------------------------------
InitUG(dim, AlgebraType("CPU", dim+1))

--------------------------------------------------------------------------------
-- Problem setup.
--------------------------------------------------------------------------------
local myProblem=require("cylinder-limex-config")
myProblem:Init(ARGS) 

--local dom = myProblem:CreateDomain(ARGS.gridName, ARGS.numRefs, ARGS.numPreRefs)
local dom = util.CreateDomain (ARGS.gridName, ARGS.numPreRefs, {"Inlet", "Outlet", "Inner", "UpperWall", "LowerWall", "CylinderWall"})
balancer.RefineAndRebalanceDomain (dom, ARGS.numRefs - ARGS.numPreRefs)

--local approxSpace = myProblem:CreateApproxSpace(dom)
-- Create the vertex-centered approximation space
approxSpace = ApproximationSpace (dom)
approxSpace:add_fct ({"u", "v", "p"}, "Lagrange", 1)

approxSpace:init_levels()
approxSpace:init_top_surface()
approxSpace:print_statistic()

util.solver.defaults.approxSpace = approxSpace

------------------------------------------------------------------------------------------
-- Compose the discretization
------------------------------------------------------------------------------------------
--local domainDisc  = myProblem:CreateDomainDisc(approxSpace)
-- inner space
NavierStokesDisc = NavierStokesFV1 ({"u", "v", "p"}, {"Inner"})
NavierStokesDisc:set_exact_jacobian (bExactJac)
NavierStokesDisc:set_stokes (bStokes)
NavierStokesDisc:set_laplace ( not(bNoLaplace) )
NavierStokesDisc:set_kinematic_viscosity (viscosity)
NavierStokesDisc:set_upwind (upwind)
NavierStokesDisc:set_peclet_blend (bPecletBlend)
NavierStokesDisc:set_stabilization (stab, diffLength)
NavierStokesDisc:set_pac_upwind (bPac)

-- boundary condition at the inlet
inflow_R = (inflow_y_1 - inflow_y_0) / 2
function inflowVel_u (x, y, t)
	return inflow * (inflow_y_1 - y) * (y - inflow_y_0) / (inflow_R * inflow_R)
end
function inflowVel2d (x, y, t)
	return inflowVel_u (x, y, t), 0.0
end
InletDisc = NavierStokesInflow (NavierStokesDisc)
InletDisc:add ("inflowVel2d", "Inlet")

-- boundary condition at the outlet
OutletDisc = NavierStokesNoNormalStressOutflow (NavierStokesDisc)
OutletDisc:add ("Outlet")

-- boundary condition at the impermeable walls
WallDisc = NavierStokesWall (NavierStokesDisc)
WallDisc:add ("UpperWall,LowerWall,CylinderWall")

-- the global discretization
domainDisc = DomainDiscretization (approxSpace)
domainDisc:add (NavierStokesDisc)
domainDisc:add (InletDisc)
domainDisc:add (OutletDisc)
domainDisc:add (WallDisc)





-- local solver = myProblem:CreateSolver(approxSpace)

--local myVelCmp = myProblem:GetVelocityCmps()
--local myDbgWriter = GridFunctionDebugWriter(approxSpace)
-- myDbgWriter.set_conn_viewer_output(true)

-- local transfer = 
-- local transfer = StdTransfer()
-- transfer:enable_p1_lagrange_optimization(true)
-- transfer:set_debug(myDbgWriter)

solverDesc =
{
	type = "newton",
	linSolver =
	{
		type = "bicgstab",
		precond =
		{
			type = "gmg",
			rap = false,
			smoother =
			{
				type = "ilu",
				beta = ilu_beta,
				consistentInterfaces = true
			},
			preSmooth = 3,
			postSmooth = 3,
			baseSolver = "lu",
			baseLevel = numPreRefs,
			approxSpace = approxSpace,
			cycle = "V"
		},
		convCheck =
		{
			type		= "standard",
			iterations	= 128,
			absolute	= 1e-8,
			reduction	= 1e-4,
			verbose		= true
		}
	},
	lineSearch =
	{
		type			= "standard",
		maxSteps		= 5,
		lambdaStart		= 1,
		lambdaReduce	= 0.5,
		acceptBest 		= true,
		checkAll		= false
	},
	convCheck =
	{
		type		= "standard",
		iterations	= 256,
		absolute	= 1e-5,
		reduction	= 1e-8,
		verbose		= true
	}
}

solver = util.solver.CreateSolver(solverDesc)
--solver:set_debug(GridFunctionDebugWriter(approxSpace))
--print(solver:config_string())


--------------------------------------------------------------------------------
-- LIMEX configuration
--------------------------------------------------------------------------------
local myLimexDesc = {}
    myLimexDesc["STAGES"] = ARGS.limexNStages -- tolerance
    myLimexDesc["TOLERANCE"] = ARGS.limexTOL-- tolerance
    myLimexDesc["REDUCTION"] = 0.5  -- reduction
    myLimexDesc["INCREASE"]  = 2.0 -- increase of time step
    myLimexDesc["SAFETY"] = 0.8    -- safety factor
    -- local errorEst = myLimexDesc["ESTIMATOR"]

    myLimexDesc["SPACES"] = {
      H1SemiComponentSpace("u", 4), 
      H1SemiComponentSpace("v", 4), 
    --  L2QuotientSpace("p", 2), -- need to factor out constants!
    }

    -- myLimexDesc["DEBUG"] = myDbgWriter -- enable output of matrices

--------------------------------------------------------------------------------
-- VTK writer.
--------------------------------------------------------------------------------
local vtkWriter = VTKOutput()
vtkWriter:clear_selection()
vtkWriter:select_all(false)
vtkWriter:select_nodal ({"u", "v"}, "vel")
vtkWriter:select_nodal ("u", "vel_u")
vtkWriter:select_nodal ("v", "vel_v")
vtkWriter:select("p", "pressure")
--vtkWriter:select(myVelCmp, "velocity")


--------------------------------------------------------------------------------
-- Solve problem.
--------------------------------------------------------------------------------
local u = GridFunction(approxSpace)
--myProblem:SetInitialValues(u,ARGS.doSteadyState)
u:set (0)
Interpolate (inflowVel_u, u, "u")
-- order the DoFs:
--OrderCuthillMcKee (approxSpace, true)
OrderLex (approxSpace,  "x")

if (ARGS.doSteadyState) then
  -- Steady state solution.
  myProblem:ComputeNonLinearSteadyStateSolution(u, domainDisc, solver)
  vtkWriter:print("CylinderSteadyState", u)
else
  -- Transient solution.
  
  if (ARGS.doInitialSteadyStateValues) then 
    myProblem:ComputeNonLinearSteadyStateSolution(u, domainDisc, solver)
  end
  
  local cTransient = 
    {
      -- Start and stop time --
      tStart = 0.0,
      tStop = 60.0,
      numTimeSteps = 60,

      -- Fractional-step-theta (same amount of work as CN, but stable) --
      --  scheme = "fracstep", maxStepSize = 0.1, minStepSize = 8.0/4096.0, redStepSize = 0.5,

      -- Limex time stepping
      scheme = "limex", maxStepSize = 0.5, minStepSize = (8.0/2048.0)*1e-8, redStepSize = 0.5, 
     
      -- Implicit Euler (requires a small time step) --
      --scheme = "impleuler", maxStepSize = 0.5, minStepSize = 0.0001, redStepSize = 0.3

      -- SDIRK (3rd-order) a.k.a. Alexander3 (we call it `sdirk3`) --
      --scheme = "sdirk", orderOrTheta = 3, maxStepSize = 0.12, minStepSize = 0.03, redStepSize = 0.5 

      -- SDIRK (4th-order) a.k.a. Hairer, Wanner, L-stable DIRK (we call it `sdirk4`), TODO: Check --
      --scheme = "sdirk", orderOrTheta = 4, maxStepSize = 0.12, minStepSize = 0.03, redStepSize = 0.5 
    }
    
    cTransient.dt = (cTransient.tStop-cTransient.tStart)/cTransient.numTimeSteps



  local ComputeEffectiveQuantities = function (u, step, time) myProblem:EvalIntegralQuantities2D(u, step, time) end

  if (cTransient.scheme == "limex") then 


      myProblem:SolveNonlinearProblemLimex(
      u, domainDisc, solver,
      vtkWriter, "limex_solution",
      cTransient.tStart, cTransient.tStop,
      cTransient.numTimeSteps, 
      cTransient.dt,  -- start step size dt
      cTransient.minStepSize, 
      cTransient.maxStepSize, 
      myLimexDesc, 
      ComputeEffectiveQuantities)
  else

    util.SolveNonlinearTimeProblem(
      u, domainDisc, solver,
      vtkWriter, "std_solution",
      cTransient.scheme, 
      cTransient.orderOrTheta, 
      cTransient.tStart, cTransient.tStop, 
      cTransient.maxStepSize, cTransient.minStepSize, cTransient.redStepSize, 
      false, false, ComputeEffectiveQuantities)

  end


end



-- Compute drag and lift coefficients.
if dim == 2 then 
  myProblem:EvalIntegralQuantities2D(u,0,0)
end	

