
local myProblem = {}


-- Reference values for Schaefer /Turek benchmarks
local ref2D_1 = {
  CD = 5.57953523384,
  CL = 0.010618948146,
  DeltaP = 0.11752016697,
  
  Um =  1.5 --0.15 -- 1.5
}

local ref2D_3 = {
  
  CD = 2.950921575, tCD=3.93625, -- maximum value and time (cited according to in John, Rang)
  CL = 0.47795, tCL= 5.693125,
  DeltaP = -0.1116,  -- at t=8
  
  Um = 1.5
}


local Cylinder2D_FE = {

  disc = {
    walls = "UpperWall,LowerWall,CylinderWall",
    inlet = function (x, y, t) return 4 * Um * y * (H-y) / (H*H), 0.0 end,
  }

}

local Cylinder3D_FE = {
 
  disc = {
    walls = "UpperWall,LowerWall,CylinderWall,FrontWall,BackWall",
    inlet = function (x, y, z, t) return 16 * Um * y * z * (H-y) * (H-z) / (H*H*H*H), 0.0, 0.0  end,
  }
  
}

GLOBAL_CYLINDER_CONFIG = { H = 0.41, L = 0.1, Um=0.5 }


function GLOBAL_CYLINDER_inletVel2d(x, y, t)
   local H  = GLOBAL_CYLINDER_CONFIG.H
   local Um = GLOBAL_CYLINDER_CONFIG.Um
   return 4 * Um * y * (H-y) / (H*H), 0.0
end

function GLOBAL_CYLINDER_inletVelX2d(x, y, t)
   local H  = GLOBAL_CYLINDER_CONFIG.H
   local Um = GLOBAL_CYLINDER_CONFIG.Um
   return 4 * Um * y * (H-y) / (H*H)
end

function GLOBAL_CYLINDER_inletVel3d(x, y, z, t)
   local H  = GLOBAL_CYLINDER_CONFIG.H
   local Um = GLOBAL_CYLINDER_CONFIG.Um
   return 16 * Um * y * z * (H-y) * (H-z) / (H*H*H*H), 0.0, 0.0
end



-- TODO: This should be integrated into a constructor!
myProblem.Init = function(self, o)

  self.dim = o.dim
  
  self.discType = o.discType
  self.vorder = o.vorder
  self.porder = o.porder
  
  
  self.doSteadyState = o.doSteadyState
  
  self.bExactJac = o.bExactJac or true
  self.bStokes= o.bStokes or false
  self.bLaplace  = o.bLaplace or false
  self.bPecletBlend  = o.bPecletBlend or false
  
  self.upwind  = o.upwind or "full"
  
  self.stab = o.stab
  self.diffLength = o.diffLength
  
  self.stabGrad = o.stabGrad
  self.stabStreamline = o.stabStreamline
  self.stabDiv =  o.stabDiv
  
  self.viscosity=1e-3  -- 1.0  -- kinematic (nu=1/Re) or dynamic (mu) does not matter, since \rho=1. 
  self.density=1.2
  self.Um = 1.5
  self.Umean2 = math.pow(2/3*self.Um, 2)
  
  self.L = GLOBAL_CYLINDER_CONFIG.L
  self.H = GLOBAL_CYLINDER_CONFIG.H
  

  
end

myProblem.CreateDomain=function(self, gridName, numRefs, numPreRefs)
  
 
  local dom = Domain()
  LoadDomain(dom, gridName)
  
  -- Create a refiner instance. This is a factory method
  -- which automatically creates a parallel refiner if required.
  local refiner =  GlobalDomainRefiner(dom)
  
  write("Pre-Refining("..numPreRefs.."): ")
  for i=1,numPreRefs do write(i .. " ");  refiner:refine(); end
  write("done. Distributing...")
  if util.DistributeDomain(dom, distributionMethod, verticalInterfaces, numTargetProcs, distributionLevel, wFct) == false then
    print("Error while Distributing Grid. Aborting.")
    exit();
  end
  write(" done. Post-Refining("..(numRefs-numPreRefs).."): ")
  for i=numPreRefs+1,numRefs do refiner:refine(); write(i-numPreRefs .. " "); end
  write("done.\n")
  
  SaveGridHierarchyTransformed(dom:grid(), dom:subset_handler(), "grid_p"..ProcRank()..".ugx", 0.5)
  
  return dom
end


myProblem.CreateApproxSpace=function (self, dom)

  local discType = self.discType
  local vorder = self.vorder
  local porder = self.porder
  
  self.approxSpace = util.ns.CreateApproxSpace(dom, discType, vorder, porder)
  
  -- print statistic on the distributed dofs
  self.approxSpace:init_levels()
  self.approxSpace:init_top_surface()
  self.approxSpace:print_statistic()
  self.approxSpace:print_local_dof_statistic(2)
  
  -- OrderLex(approxSpace,"x")
  -- OrderCuthillMcKee(approxSpace,true)
  return self.approxSpace
end

-- Extract velocity components from approx space.
myProblem.GetVelocityCmps =function(self)
local FctCmp = self.approxSpace:names()
local VelCmp = {}
for d = 1,#FctCmp-1 do VelCmp[d] = FctCmp[d] end
 return VelCmp
end

--------------------------------------------------------------------------------
-- Discretization
--------------------------------------------------------------------------------


-- Creates the domain discretization

myProblem.CreateDomainDisc = function (self,approxSpace) --, vorder, porder)
  local discType = self.discType
  local FctCmp = approxSpace:names()
  local NavierStokesDisc = NavierStokes(FctCmp, {"Inner"}, discType)
  NavierStokesDisc:set_exact_jacobian(self.bExactJac)
  NavierStokesDisc:set_stokes(self.bStokes)
  NavierStokesDisc:set_laplace(not(self.bLaplace))
  NavierStokesDisc:set_kinematic_viscosity( self.viscosity );
  --NavierStokesDisc:set_density( self.density );
        
  local dim = self.dim      
  local vorder = approxSpace:lfeid(0):order()
  local porder = approxSpace:lfeid(self.dim):order()
 
  
  --upwind if available
  if discType == "fv1" or discType == "fvcr" then
    NavierStokesDisc:set_upwind(self.upwind)
    NavierStokesDisc:set_peclet_blend(self.bPecletBlend)
  end
  
  -- fv1 must be stablilized
  if discType == "fv1" then
    NavierStokesDisc:set_stabilization(self.stab, self.diffLength)
    NavierStokesDisc:set_pac_upwind(false)
  end
  
  -- fe must be stabilized for (Pk, Pk) space
  if (discType == "fe") and (porder == vorder) then
    NavierStokesDisc:set_stabilization(self.stabGrad)
    NavierStokesDisc:set_stab_streamline(self.stabStreamline) 
  end
  
  if discType == "fe" then
   NavierStokesDisc:set_quad_order(math.pow(vorder, dim)+2)
   -- NavierStokesDisc:set_quad_order(3)
   NavierStokesDisc:set_stab_div(self.stabDiv)
  end
  if discType == "fv" then
    NavierStokesDisc:set_quad_order(math.pow(vorder, dim)+2)
  end
  
  -- setup Outlet
   OutletDisc = NavierStokesNoNormalStressOutflow(NavierStokesDisc)
   OutletDisc:add("Outlet")
  
  -- setup Inlet

  
  local InletDisc = NavierStokesInflow(NavierStokesDisc)
  InletDisc:add("GLOBAL_CYLINDER_inletVel"..dim.."d", "Inlet")
  
  -- John's (physically unrealistic) BC
  --local InletDisc = DirichletBoundary()
  --InletDisc:add("GLOBAL_CYLINDER_inletVelX"..dim.."d", "u" ,"Inlet, Outlet")
  --InletDisc:add(0.0, "v" ,"Inlet, Outlet")
  
  --setup Walls
  local WallDisc = NavierStokesWall(NavierStokesDisc)
  if dim == 2 then
    WallDisc:add("UpperWall,LowerWall,CylinderWall")
  elseif dim == 3 then
    WallDisc:add("UpperWall,LowerWall,CylinderWall,FrontWall,BackWall") 
  end

  local DirichletBnd = DirichletBoundary()
  -- DirichletBnd:add(0, "p", "FIXP")  -- fix pressure => oscillations
  -- DirichletBnd:add("inletVelX2d", "u", "Outlet")  -- fix pressure
  -- DirichletBnd:add(0.0, "v", "Outlet")  -- fix pressure
  
  --local NoSlipBnd = DirichletBoundary()
  --NoSlipBnd:add(0, "u", "UpperWall,LowerWall,CylinderWall")  -- no slip
  --NoSlipBnd:add(0, "v", "UpperWall,LowerWall,CylinderWall")  -- no slip
  
  
  -- Finally we create the discretization object which combines all the
  -- separate discretizations into one domain discretization.
  local domainDisc = DomainDiscretization(approxSpace)
  domainDisc:add(NavierStokesDisc)
  domainDisc:add(InletDisc)
  domainDisc:add(WallDisc)
  --domainDisc:add(DirichletBnd)
  --domainDisc:add(NoSlipBnd)
  domainDisc:add(OutletDisc)
  
  return domainDisc
end

--------------------------------------------------------------------------------
-- Solution of the Problem
--------------------------------------------------------------------------------
myProblem.CreateSolver = function (self, approxSpace, discType, p)

  local discType=self.discType
  local p = nil
  
  local base = LU()
  
  local smoother = nil
  if discType == "fvcr" or discType == "fecr" then 
    smoother = ComponentGaussSeidel(0.1, {"p"}, {1,2}, {1})
  elseif discType == "fv1" then 
    smoother = ILU()
    smoother:set_damp(0.7)
  else
     smoother = ComponentGaussSeidel(1.0, {"p"})
     smoother:set_alpha(1.0)
     smoother:set_beta(1.0)
     smoother:set_weights(true)
  end
  
  local numPreSmooth, numPostSmooth, baseLev, cycle, bRAP = util.gmg.parseParams()
  local cycleType = "W"
  local bRAP= true
  local gmg = util.gmg.create(approxSpace, smoother, numPreSmooth, numPostSmooth,
               cycleType, base, baseLev, bRAP)
  --gmg:add_prolongation_post_process(AverageComponent("p"))
  local transfer = StdTransfer()
  transfer:enable_p1_lagrange_optimization(false)
  gmg:set_transfer(transfer)
  
  local sol = util.solver.parseParams()
  local solver = util.solver.create(sol, gmg)
  if self.bStokes then
    solver:set_convergence_check(ConvCheck(10000, 5e-12, 1e-99, true))
  else 
    solver:set_convergence_check(ConvCheck(10000, 5e-12, 1e-2, true)) 
  end
    
  local convCheck = ConvCheck(50, 1e-11, 1e-99, true)
  
  local newtonSolver = NewtonSolver()
  newtonSolver:set_linear_solver(solver)
  newtonSolver:set_convergence_check(convCheck)
  newtonSolver:set_line_search(StandardLineSearch(10, 1.0, 0.9, true, true))
  newtonSolver:set_debug(GridFunctionDebugWriter(approxSpace))
  
  return newtonSolver
end


myProblem.ComputeNonLinearSteadyStateSolution = function(self, u, domainDisc, solver, cmp)
  --local cmp = cmp or "p"
  --util.rates.static.StdComputeNonLinearSolution(u, domainDisc, solver)
  --AdjustMeanValue(u, cmp)
  
  solver:init(AssembledOperator(domainDisc))
  solver:prepare(u)

  -- apply the solver for the stationary pressure problem
  tBefore_s= os.clock()
  if not solver:apply(u) then
  	print("===> THE PREPARATION PHASE FAILED! <===")
  	exit()
  end
	
  tAfter_s = os.clock()
  print("Computation for steady state took " .. tAfter_s-tBefore_s .. " seconds.")

end

myProblem.SetInitialValues = function (self, u)
  u:set(0.0)
end

-- Evaluate drag, lift and deltaP
myProblem.EvalIntegralQuantities2D = function (self, u, step, time)
    print("EvalIntegralQuantities2D")

    local DL = DragLift(u, "u,v,p", "CylinderWall", "Inner", self.viscosity, 1.0, self.vorder+3)
    
    local L = self.L
    
    local C_D = 2*DL[1]/(self.Umean2*L)
    local C_L = 2*DL[2]/(self.Umean2*L)
  
    local PEval = GlobalGridFunctionNumberData(u, "p")
    local Delta_P = PEval:evaluate(Vec2d(0.15, 0.2)) - PEval:evaluate(Vec2d(0.25, 0.2))
  
    print("EVAL_P1:\t"..time.."\t"..PEval:evaluate(Vec2d(0.15, 0.2)))
    print("EVAL_P2:\t"..time.."\t"..PEval:evaluate(Vec2d(0.25, 0.2)))
    print("EVAL_DELTA_P:\t"..time.."\t"..Delta_P)
  
    print("EVAL_C_D:\t"..time.."\t"..C_D)
    print("EVAL_C_L:\t"..time.."\t"..C_L)
end

myProblem.SolveNonlinearProblemLimex = function (
  self,
  u,
  domainDisc,
  newtonSolver,
  out,
  filename,
  startTime,
  endTime,
  numTimeSteps,
  dt,
  minStepSize,
  maxStepSize,
  adaptiveDesc,
  postProcess)


  -- read adaptive stuff
  local inc_fac = adaptiveDesc["INCREASE"] or 1.5        -- increase of time step
 
  -- check parameters
  if filename == nil then filename = "sol" end
  if minStepSize == nil then minStepSize = maxStepSize end


  -- Check input parameters.
  if u == nil or domainDisc == nil or newtonSolver == nil
    or startTime == nil or endTime == nil or maxStepSize == nil then
    print("Wrong usage found. Please specify parameters as below:")
    
    if (u == nil) then print ("Did not find u!"); end;
    if (domainDisc == nil) then print ("Did not find domainDisc!"); end;

    if (startTime == nil) then print ("Did not find endTime!"); end;
    if (endTime == nil) then print ("Did not find endTime!"); end;
    if (maxStepSize == nil) then print ("Did not find maxStepSize!"); end;
    --util.PrintUsageOfSolveTimeProblem()
    exit()
  end


  print ("maxStepSize ="..maxStepSize)
  print ("minStepSize ="..minStepSize)

  print ("startTime ="..startTime)
  print ("endTime ="..endTime)
  
  
  -- Create LIMEX descriptor
  local limexDesc = {

        nstages = adaptiveDesc["STAGES"] or 2,
        steps = {1,2,3,4,5,6},
        nthreads = 1, 
        tol = adaptiveDesc["TOLERANCE"] or  1e-3,
        rhoSafetyOPT = adaptiveDesc["SAFETY"] or 0.25,

        dt = dt,
        dtmin = minStepSize,
        dtmax = maxStepSize,
        dtred = adaptiveDesc["REDUCTION"] or 0.5,  -- reduction of time step

        -- set disc & solver
        domainDisc= domainDisc,
        nonlinSolver = newtonSolver,   
        -- makeConsistent = true,

        matrixCache = true, -- or true,
        -- costStrategyOPT = time.limexDesc.costStrategyOPT,
        debugOPT = 5,

       -- dampScheideggerOPT = time.limexDesc.dampScheideggerOPT or 1.0,
       -- partialVeloMaskOPT = time.limexDesc.partialVeloMaskOPT or 0,
     }
  -- Create LIMEX object
  local limex = util.limex.CreateIntegrator(limexDesc)
      
   limex:set_time_step(limexDesc.dt)
   limex:set_dt_min(limexDesc.dtmin)
   limex:set_dt_max(limexDesc.dtmax)
   limex:set_reduction_factor(limexDesc.dtred)
  
    if (adaptiveDesc["DEBUG"]) then 
      --limex:set_debug(adaptiveDesc["DEBUG"])
      limex:set_debug_for_timestepper(adaptiveDesc["DEBUG"])
    end
  -- Register LUA callback.
  if type(postProcess) == "function" then 
    -- a) LUA functions
    local luaobserver = LuaCallbackObserver()
    
     function __util_LimexLuaCallbackPost(step, t, currdt) 
            local sol=luaobserver:get_current_solution()
            print(postProcess)
            postProcess(sol, step, t, currdt)
            return 1
      end
      luaobserver:set_callback("__util_LimexLuaCallbackPost") 
    limex:attach_observer(luaobserver)
   end
   
   
   local limexErrorEst 
   limexErrorEst = CompositeGridFunctionEstimator() 
   
   if (type(adaptiveDesc["SPACES"])=="table") then
    for i, _spacei in ipairs(adaptiveDesc["SPACES"]) do 
      print(_spacei)
      limexErrorEst:add(_spacei)
      
    end
   end -- table
  
 
  -- limex:set_space(limexErrorSpace)


   print(limexErrorEst:config_string())
   limex:add_error_estimator(limexErrorEst)
   
   -- Solve problem
    print(">> Solve using LIMEX...")
    
    -- Replace convergence check.
    local limexConvCheck = ConvCheck()
    limexConvCheck:set_maximum_steps(10)
    limexConvCheck:set_minimum_defect(1e-5)
    limexConvCheck:set_reduction(1e-9)
    limexConvCheck:set_verbose(true)
    --limexConvCheck:set_supress_unsuccessful(true)
    
    newtonSolver:set_convergence_check(limexConvCheck) 
    newtonSolver:disable_line_search()
    
    -- Execute solver
    local sw = CuckooClock()
    sw:tic()
    out:print("Sol", u,0,startTime)	
    print(newtonSolver:config_string())
    for step = 1, numTimeSteps do
	StepTime=startTime+(step)*dt
	limex:apply(u, StepTime, u, StepTime-dt)
	out:print("Sol", u,step,StepTime)	
    end	
    print ("CDELTA="..sw:toc())
  return 
end

return myProblem
    