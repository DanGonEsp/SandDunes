--------------------------------------------------------------------------------
-- Initial Parameters
--------------------------------------------------------------------------------
local myProblem = {}

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
	
	
	self.boolSolverDesc = false
	self.NewtonSolverDescSteady = {}
	self.NewtonSolverDesc = {}
  
  
end



--------------------------------------------------------------------------------
-- SOLVER
--------------------------------------------------------------------------------
myProblem.CreateSolver = function (self, domainDisc, approxSpace, timeDisc)

		-- For debugging only (to write the intermediate data): --
	--util.debug_dir = "FLOW_DEBUG"
	--util.debug = { vtk = true, conn_viewer = false }
	--util.CreateGridFuncDebugWriter (approxSpace)
	-- --
	
	----------------------------------------------------------
	-- preconditioners used by FETI sub solvers
	----------------------------------------------------------
	-- base solver


	basePre = ILUT()
	basePre:set_threshold(1e-7)

	--basePre=ElementGaussSeidel()
	
	
	baseSolver = LinearSolver()
	baseSolver:set_preconditioner(basePre)
	baseSolver:set_convergence_check(ConvCheck(10000, 1e-07,1e-04,true))
	
	--baseSolver=LU()
	
	
	
	local LinearSolverDesc =
	{
		type = "bicgstab",
		precond =
		{
			type = "gmg",
			rap = true,
			rim = false,
			cycle = "V",
			smoother =
			{
				type = "ilu",
				beta = self.value_beta,
				damping 	= self.damping_mg,
				sort	= false,
				--sortEps 	= 1.e-50,
				inversionEps 	= 1.e-16,
				consistentInterfaces   = false,     --consistentInterfaces and overlap shouldnot be activated at the same time
				overlap 		= true,             --consistentInterfaces and overlap shouldnot be activated at the same time
				--ordering 		= nil
			},
			preSmooth = 2,
			postSmooth = 2,
			baseSolver = ilutSolver,
			baseLevel = self.numPreRefs
		},
		convCheck =
		{
			type		= "standard",
			iterations	= self.max_linear_steps,
			absolute	= 1e-7,
			reduction	= 1e-3,
			verbose		= true
		}
	}
	NewtonConvCheckSteady =
	{
		type		= "standard",
		iterations	= self.max_newton_steps_steady_state,
		absolute	= self.AbsDefect,
		reduction	= self.RedDefect,
		verbose		= true
	}
	NewtonConvCheck =
	{
		type		= "standard",
		iterations	= self.max_newton_steps_transcient,
		absolute	= self.AbsDefect,
		reduction	= self.RedDefect,
		verbose		= true
	}
	NewtonSolverDescSteady =
	{
		type = "newton",
		debug = false, -- for the debug output from the Newton's method
		linSolver = LinearSolverDesc
		,
		lineSearch =
		{
			type			= "standard",
			maxSteps		=self.lambdamaxSteps,
			lambdaStart		= self.lambdaStart,
			lambdaReduce	= 0.7,
			acceptBest 		= true,
			checkAll		= false,
			suffDesc		= 0.3,
			maxDefect	= 2e20
			
		},
		convCheck = NewtonConvCheckSteady
	}

	NewtonSolverDesc = NewtonSolverDescSteady
	NewtonSolverDesc.convCheck = NewtonConvCheck

	local NewtonSolverSteady = nil
	if self.doSteadyState then NewtonSolverSteady = util.solver.CreateSolver(NewtonSolverDesc) end
	local LinearSolver = util.solver.CreateSolver(LinearSolverDesc)

	local limexConvCheck=ConvCheck(1, 5e-5, 1e-8, true)
	limexConvCheck:set_supress_unsuccessful(true)



	local NLSolver = nil
	local limex = nil
	if self.timeMethod == "limex" then
		NLSolver = NewtonSolver()
		NLSolver:set_linear_solver(LinearSolver)
		NLSolver:set_convergence_check(limexConvCheck)
		limex = myProblem:LimexObject( domainDisc, NLSolver)
	else
		NLSolver = util.solver.CreateSolver(NewtonSolverDesc)
		op = AssembledOperator(timeDisc)
		op:init()
		NLSolver:init(op)
		if NLSolver:prepare(u) == false then
			print ("Newton solver prepare failed.") exit()
		end
	end

	self.NewtonSolverDescSteady = NewtonSolverDescSteady
	self.NewtonSolverDesc = NewtonSolverDesc
	self.boolSolverDesc = true

	  
	return op, NLSolver, NewtonSolverSteady, limex
end



--------------------------------------------------------------------------------
-- Writing Output parameters
--------------------------------------------------------------------------------

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


--------------------------------------------------------------------------------
-- SteadyState Solution
--------------------------------------------------------------------------------
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

--------------------------------------------------------------------------------
-- Solution of NonLinear Problem  (Euler temporal discretization)
--------------------------------------------------------------------------------
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

	if(self.boolSolverDesc == false) then print("SolverDesc Not initialized ") exit() end
	local solverDesc = self.NewtonSolverDesc
		
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
							absolute	= self.AbsDefect,
							reduction	= self.RedDefect,
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

--------------------------------------------------------------------------------
-- LIMEX
--------------------------------------------------------------------------------
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

	local dtmax = self.DTmax
	local dtmin = self.DTmin
	local dtlimex = 1.0
	local gridSize = 1.0
	--  Euclidean (algebraic) norm
	--local estimator = Norm2Estimator()
	--tol = 0.37/(gridSize)*tol


	--print (estimator)
	local limexEstimator = CompositeGridFunctionEstimator()
	--limexEstimator:add(H1SemiComponentSpace("u", 2 ))
	--limexEstimator:add(H1SemiComponentSpace("v", 2 ))

	limexEstimator:add(L2ComponentSpace("u", 2))
	limexEstimator:add(L2ComponentSpace("v", 2))

	limexEstimator:add(H1SemiComponentSpace("p", 2))
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
	  dtmax = dtmax,
	  dtmin = dtmin,
	  rhoSafetyOPT = 0.2,
	  dtred = 0.5,
	  dtIncr = 1.3,
	  matrixCache = true,
	  conservative = true
	  
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
	print ("tolLimex  = "..self.tol)
	return limex
end
--------------------------------------------------------------------------------
-- SolutionNonLinearProblem LIMEX
--------------------------------------------------------------------------------
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
    
