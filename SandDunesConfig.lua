
local myProblem = {}

-- TODO: This should be integrated into a constructor!
myProblem.Init = function(self, o)

		-- Numerical parameters of the discretization
	self.dim = o.dim
	self.file_name = o.file_name
	self.elem_type = o.elem_type
	self.numRefs = o.numRefs
	self.numPreRefs = o.numPreRefs
	self.startTime = o.startTime
	self.endTime = o.endTime
	self.numTimeSteps = o.numTimeSteps
	self.DTmax = o.DTmax
	self.DTmin = o.DTmin
	self.outputFactor = o.outputFactor

	self.timeMethod  = o.timeMethod
	self.modifyDT  = o.modifyDT
	self.incr_factor = o.incr_factor
	self.red_factor_fail = o.red_factor_fail
	self.red_factor_success = o.red_factor_success
	self.optimal_newton_steps = o.optimal_newton_steps
	self.maxConvRate = o.maxConvRate
	self.minConvRate = o.minConvRate
	self.boolDebugStep = o.boolDebugStep
	
	self.tol = o.tol
	self.nstages = o.nstages
	self.limex_partial_mask = o.limex_partial_mask
	self.limex_debug_level = o.limex_debug_level

	self.max_newton_steps_steady_state = o.max_newton_steps_steady_state
	self.max_newton_steps_transcient = o.max_newton_steps_transcient
	self.AbsDefect = o.AbsDefect
	self.RedDefect = o.RedDefect

	self.max_linear_steps = o.max_linear_steps
	self.damping_mg = o.damping_mg
	self.value_beta = o.value_beta

	self.lambdamaxSteps = o.lambdamaxSteps
	self.lambdaStart = o.lambdaStart
	
		-- Physical phenomenon of simulation
	self.doSteadyState = o.doSteadyState
	self.boolSource = o.boolSource
	self.consistentRho_in_source = o.consistentRho_in_source
	self.boolRelativeVel = o.boolRelativeVel
	self.boolGradientPsSource = o.boolGradientPsSource
	self.boolViscPs = o.boolViscPs
	self.boolAveDiff = o.boolAveDiff
	self.boolSlipDiff = o.boolSlipDiff
	
	self.inflow = o.inflow
	self.H_0 = o.H_0
	self.ReferencePressure = o.ReferencePressure
	self.bStokes = o.bStokes or false
	self.bNoLaplace = o.bNoLaplace or false
	self.bExactJac = o.bExactJac or false
	self.bPecletBlend = o.bPecletBlend or false
	self.upwind_m = o.upwind_m or "full"
	self.upwind_t = o.upwind_t or "full"
	self.upwind_r = o.upwind_r or "full"
	self.bPac = o.bPac
	self.diffLength = o.diffLength
	self.stab = o.stab
	self.turbViscMethod = o.turbViscMethod
	self.modellconstant = o.modellconstant

	
	--Material Properties
	self.nu_a = o.nu_a
	self.rho_a = o.rho_a
	self.rho_s = o.rho_s
	self.dp = o.dp
	self.nu_s = o.nu_s
	self.c_init = o.c_init

	self.alpha_max = o.alpha_max
	self.alpha_min = o.alpha_min
	self.granular_model = o.granular_model
	self.density_model = o.density_model
	self.interface_value = o.interface_value
	self.drag_mod = o.drag_mod


	self.FR = o.FR
	self.B_phi = o.B_phi
	self.deltaGamma = o.deltaGamma
	self.Visc_limit = o.Visc_limit

	self.deltaPs = o.deltaPs
	self.deltaI = o.deltaI
	self.FricMu_1 = o.FricMu_1
	self.FricMu_2 = o.FricMu_2
	self.I_0 = o.I_0
  
  
end


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

myProblem.WriteValues = function (self, folder, step, time, Value_inner1, Value_inner2, WorkTime, Newton_Steps, Newton_Steps_fail,linsolver_calls,linsolver_steps,boolTotal)
	if(boolTotal) then
		file = io.open(folder .. "/Integral.txt", "a")
		file:write(string.format("-----------------------------------------------------------------------------------------------------------\n"))
		file:close()
	end
	if(step == 0) then
		file = io.open(folder .. "/Integral.txt", "w+")
		file:write(string.format("Step\tTime\t\tVol-Dom_1\tVol-Dom_2\tWork time\tTNSteps\tSNSteps\tFNSteps\tLinCalls LinSteps\n"))
		file:write(string.format("%d\t%.6f\t%.6f\t%.6f\t%.6f\t%d\t%d\t%d\t%d\t%d\n", step, time, Value_inner1, Value_inner2, WorkTime, Newton_Steps, Newton_Steps-Newton_Steps_fail, Newton_Steps_fail,linsolver_calls,linsolver_steps))
		file:close()
	else
		file = io.open(folder .. "/Integral.txt", "a")
		file:write(string.format("%d\t%.6f\t%.6f\t%.6f\t%.6f\t%d\t%d\t%d\t%d\t%d\n", step, time, Value_inner1, Value_inner2, WorkTime, Newton_Steps, Newton_Steps-Newton_Steps_fail, Newton_Steps_fail,linsolver_calls,linsolver_steps))
		file:close()
	end
end

myProblem.ComputeNonLinearSteadyStateSolution = function(self, u, domainDisc, solver)

	-- Fix the mass fraction and solve the linear problem for the momentum
	fixer = DirichletBoundary()
	domainDisc:add(fixer)
	fixer:invert_subset_selection()
	fixer:add("c", "")

	solver:init(AssembledOperator(domainDisc))
	
	solver:prepare(u)

	-- apply the solver for the stationary pressure problem
     print("++++++ STEADY STATE CALCULATION BEGIN ++++++")
	tBefore_s= os.clock()
	if not solver:apply(u) then
		print("===> THE PREPARATION PHASE FAILED! <===")
		exit()
	end
	tAfter_s = os.clock()
    num_newton_steps = solver:num_newton_steps()
    linsolver_calls = solver:total_linsolver_calls()
    linsolver_steps = solver:total_linsolver_steps()
    average_linear_steps = solver:total_average_linear_steps()
    average_non_linear_rates = solver:total_average_non_linear_rates()
    
    print("num_newton_steps = " .. num_newton_steps .. ".")
    print("linsolver_calls = " .. linsolver_calls .. ".")
    print("linsolver_steps = " .. linsolver_steps .. ".")
    print("average_linear_steps = " .. average_linear_steps .. ".")
    print("average_non_linear_rates = " .. average_non_linear_rates .. ".")
    
    solver:clear_average_convergence();
    
	time_work_steady = tAfter_s-tBefore_s
	print("Computation for steady state took " .. time_work_steady .. " seconds.")
	domainDisc:remove (fixer)
	print("++++++++++++++++++++++++ INITIAL CONDITIONS  (STEADY STATE DONE) ++++++++++++++++++++++++")
	return time_work_steady, linsolver_calls, linsolver_steps

end


myProblem.SolveNonlinearProblem = function (self, u, solver, op, timeDisc, solTimeSeries, dt, step, StartTime, EndTime)

	                    
	local red_factor_fail = self.red_factor_fail
	local red_factor_success = self.red_factor_success
	local incr_factor = self.incr_factor
	local optimal_newton_steps = self.optimal_newton_steps
	local DTmin = self.DTmin
	local DTmax = EndTime-StartTime
	local modifyDT = self.modifyDT
	local boolDebugStep = self.boolDebugStep
	local AbsDefect = self.AbsDefect
	local RedDefect = self.RedDefect
	local maxConvRate = self.maxConvRate
	local minConvRate = self.minConvRate

	local solverDesc = self.solverDesc
		
    local Newton_Steps = 0
    local Newton_Steps_fail = 0
    local linsolver_calls_step = 0
    local linsolver_steps_step = 0
    local CompletedStep = false
    local time = StartTime
    local time2 = StartTime
    local dt_in = dt

        
	while CompletedStep==false  do
        Newton_Steps = Newton_Steps+1
		-- choose time step

		do_dt = math.min(dt_in,math.max((time+DTmax-time2), 0.0))
		print("Size of timestep dt: " .. do_dt)
		-- setup time Disc for old solutions and timestep
		timeDisc:prepare_step(solTimeSeries, do_dt)
	
		-- prepare newton solver
		if solver:prepare(u) == false then
			print ("Newton solver failed at step "..step.."."); exit();
		end
	
		-- apply newton solver
            
        print("++++++ TIMESTEP " .. step-1 + (time2+do_dt-time)/DTmax .. " BEGIN ++++++")
		if solver:apply(u)  == false then
            Newton_Steps_fail = Newton_Steps_fail+1
            dt_in = math.max(dt_in*red_factor_fail,0.99999*DTmin)
            print("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<")
            print("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<   Reducing Timestep in step  " .. step-1 + (time2+do_dt-time)/DTmax .. ".")
            print("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<   new DT          =    " .. dt_in .. "     Time = " .. time2 .."")
			print("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<   old DT    		=    " .. do_dt .."")
            print("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<   DTmax       =    " .. DTmax .. "     DTmin = " .. DTmin .."")
            print("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<   Rho         =    " .. solver:total_average_non_linear_rates() .. ".")
            print("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<   Dt_factor   =    " .. red_factor_fail .. ".")
            
			
			if dt_in < DTmin  or modifyDT== false then
				print ("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx            Time step below minimum. Aborting. Failed at step  " .. step-1+ (time2+do_dt-time)/DTmax .. ".");
				if boolDebugStep then
					convCheck =
						{
							type		= "standard",
							iterations	= 5,
							absolute	= AbsDefect,
							reduction	= RedDefect,
							verbose		= true
						}

					solverDesc.debug = true
					solverDesc.convCheck = convCheck
					solver = util.solver.CreateSolver(solverDesc)
					solver:init(op)
					timeDisc:prepare_step(solTimeSeries, do_dt)
					if solver:prepare(u) == false then
						print ("Newton solver failed at DEBUG step "..step..".");
					end
					print("++++++ DEBUG STEP  BEGIN ++++++")
					solver:apply(u)
					print ("Newton solver failed at DEBUG step "..step..".");
				end
				
				exit();

			else
				VecScaleAssign(u, 1.0, solTimeSeries:latest())
			
			end
		else
            
            if  (time2 + do_dt + (1e-07) * DTmin -time)/DTmax > 1.0 then
            
                time= timeDisc:future_time()                                        -- update new time
                                                                                
                oldestSol = solTimeSeries:oldest()                                  -- get oldest solution
                                                                                    
				VecScaleAssign(oldestSol, 1.0, u)                                             -- copy values into oldest solution (we reuse the memory here)
                                                                                    
                solTimeSeries:push_discard_oldest(oldestSol, time)                  -- push oldest solutions with new values to front, oldest sol pointer is poped from end
                
                CompletedStep = true
                
                
            else
                
                time2 = timeDisc:future_time()                                      -- update new time
                                                                                    
                oldestSol = solTimeSeries:oldest()                                  -- get oldest solution
                                                                                    
				VecScaleAssign(oldestSol, 1.0, u)                                             -- copy values into oldest solution (we reuse the memory here)
                                                                                    
                solTimeSeries:push_discard_oldest(oldestSol, time2)                 -- push oldest solutions with new values to front, oldest sol pointer is poped from end
                
                CompletedStep = false
                
            
            end
            
            average_non_linear_rates = solver:total_average_non_linear_rates()
            num_newton_steps = solver:num_newton_steps()
            
            if CompletedStep then
                frac_step = 0
                Local_Time = time
                Dt_factor = dt_in/do_dt
            else
                frac_step = -1 + (time2-time)/DTmax
                Local_Time = time2
            end
            
            
            if modifyDT then
                if(average_non_linear_rates>maxConvRate and num_newton_steps>30) then
                    dt_in = math.max(do_dt*red_factor_success,1.00001*DTmin)
                    print ("-------------------------------------------------------------------Time step decrease at Step " .. step-1 + (time2-time)/DTmax .. ", dt =  " .. dt_in .. ". ")
                else if (CompletedStep== false or Dt_factor < 0.98 ) then
                    if(average_non_linear_rates<minConvRate or num_newton_steps<optimal_newton_steps) then
                        dt_in=math.min(incr_factor*dt_in,DTmax)
                        print ("-------------------------------------------------------------------Time step increased at Step " .. step-1 + (time2-time)/DTmax ..", dt =  " .. dt_in .. ". ")
                    end
                    end
                end
            end
            
            
            
            CFL=cflNumber(u,do_dt)                                              -- compute CFL number
            print("DT=" .. do_dt .. "");
            print("Time=" .. Local_Time .. "");
            
            print("++++++ TIMESTEP " .. step + frac_step.. "  END ++++++")
            print("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<")
            print("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<   Successful semi Timestep in step  " .. step + frac_step  .. ".")
            print("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<   new DT    =    " .. dt_in .. "     Time = " .. Local_Time .."")
            print("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<   old DT    =    " .. do_dt .."")
            print("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<   DTmax =    " .. DTmax .. "     DTmin = " .. DTmin .."")
            print("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<   Dt_factor   =    " .. dt_in/do_dt .. ".")
            print("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<   Rho   =    " .. average_non_linear_rates .. ".")
                        

		end
  
        
        linsolver_calls_step = linsolver_calls_step + solver:total_linsolver_calls()
        linsolver_steps_step = linsolver_steps_step + solver:total_linsolver_steps()
        solver:clear_average_convergence();
			
	end

  return Newton_Steps, Newton_Steps_fail, linsolver_calls_step, linsolver_steps_step, dt_in
end

myProblem.LimexObject = function ( self, domainDisc, limexNLSolver)


	-- local refObserver = PlotRefOutputObserver("DirichletValue", vtk) -- now obsolete
	local luaObserver = LuaCallbackObserver()

	-- work-around (waiting for implementation of SmartPtr forward to lua...)
	function luaPostProcess(step, time, currdt)
	  print("LUAPostProcess: "..step..","..time..","..currdt)
	  postProcess(luaObserver:get_current_solution(), step, time, currdt)
	  return 0;
	end
	luaObserver:set_callback("luaPostProcess")

	local dtlimex = self.DTmax
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

	  nstages = self.nstages,
	  steps = {1,2,3,4,5,6},
	  nthreads = 1,
	  domainDisc=domainDisc,
	  nonlinSolver = limexNLSolver,
	  -- makeConsistent = true,
	  
	  tol = self.tol,
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
	print ("tolLimex  = "..params.tol)
	return limex
end

myProblem.SolveNonlinearProblemLimex = function (self, u, limex, NLSolver, StartTime, EndTime)

	limex:apply(u, EndTime, u, StartTime)
	
    local Newton_Steps = NLSolver:total_linsolver_calls()
	local Newton_Steps_fail = 0
	linsolver_calls_step = NLSolver:total_linsolver_calls()
	linsolver_steps_step =  NLSolver:total_linsolver_steps()
	NLSolver:clear_average_convergence();

  return Newton_Steps, Newton_Steps_fail, linsolver_calls_step, linsolver_steps_step
end

return myProblem
    
