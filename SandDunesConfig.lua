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
	self.DT = o.DT
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
	self.NewtonDebug = o.NewtonDebug
	self.NewtonSteadyDebug = o.NewtonSteadyDebug
	self.NewtonUpdater = o.NewtonUpdater
	self.StepDebug = o.StepDebug
	self.debug_dir = o.debug_dir
	
		
	self.tol = o.tol
	self.nstages = o.nstages
	self.limex_partial_mask = o.limex_partial_mask
	self.limex_debug_level = o.limex_debug_level

	self.max_newton_steps_steady_state = o.max_newton_steps_steady_state
	self.max_newton_steps_transient = o.max_newton_steps_transient
	self.AbsDefect = o.AbsDefect
	self.RedDefect = o.RedDefect

	self.LinAbsDefect = o.LinAbsDefect
	self.LinRedDefect = o.LinRedDefect
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
	
	----------------------------------------------------------
	-- LineSearch
	----------------------------------------------------------
	local NewtonLineSearch = nil
	if true then
		NewtonLineSearch = StandardLineSearch()
		NewtonLineSearch:set_maximum_steps(self.lambdamaxSteps)
		NewtonLineSearch:set_lambda_start(self.lambdaStart)
		NewtonLineSearch:set_reduce_factor(0.5)
		NewtonLineSearch:set_accept_best(true)
		NewtonLineSearch:set_check_all(false)
		NewtonLineSearch:set_suff_descent_factor(0.2)
		NewtonLineSearch:set_maximum_defect(2e20)
	else
		NewtonLineSearch = TrustRegionMethod()
		NewtonLineSearch:set_maximum_steps(3)
		NewtonLineSearch:set_lambda_start(self.lambdaStart)
		NewtonLineSearch:set_reduce_factor(0.5)
		NewtonLineSearch:set_accept_best(true)
		NewtonLineSearch:set_check_all(false)
		NewtonLineSearch:set_suff_descent_factor(0.25)
		NewtonLineSearch:set_maximum_defect(2e20)
	end

	----------------------------------------------------------
	-- NoLinear COnvCheck
	----------------------------------------------------------
	local NewtonSteadyConvCheck=ConvCheck(self.max_newton_steps_steady_state, self.AbsDefect, self.RedDefect, true)
	local NewtonConvCheck=ConvCheck(self.max_newton_steps_transient, self.AbsDefect, self.RedDefect, true)
	local LinearConvCheck=ConvCheck(self.max_linear_steps, self.LinAbsDefect, self.LinRedDefect, true)
	local LimexConvCheck=ConvCheck(1, self.AbsDefect, 1e-8, true)
	      LimexConvCheck:set_supress_unsuccessful(true)
	
	----------------------------------------------------------
	-- Smoothers
	----------------------------------------------------------
	-- base solver
	baseSolver = LU()
	baseSolver = AgglomeratingSolver(SuperLU());
	
	
	ilu = ILU()
	ilu:set_beta(self.value_beta)
	ilu:set_damp(self.damping_mg)
	ilu:set_ordering_algorithm(TopologicalOrdering())
	ilu:set_sort(true)
	--ilu:set_sort_eps(1.e-50)
	ilu:set_inversion_eps(1.e-16)
	ilu:enable_consistent_interfaces(true)
	ilu:enable_overlap(false)
	
	jac = Jacobi (0.7);
	
	bgs = BlockGaussSeidel ();
	
	gs = GaussSeidel()
	gs:enable_consistent_interfaces(false)
	gs:enable_overlap(false)
	
	sgs = SymmetricGaussSeidel ()
	sgs:enable_consistent_interfaces(true)
	sgs:enable_overlap(false)

	egs = ElementGaussSeidel();

	cgs = ComponentGaussSeidel(0.1, {"p"}, {1,2}, {1})
	

	----------------------------------------------------------
	-- preconditioners
	----------------------------------------------------------

	gmg = GeometricMultiGrid(approxSpace)
	gmg:set_discretization(domainDisc)
	gmg:set_base_level(self.numPreRefs)
	gmg:set_base_solver(baseSolver)
	gmg:set_smoother(ilu)
	gmg:set_cycle_type(1)
	gmg:set_num_presmooth(2)
	gmg:set_num_postsmooth(2)
	gmg:set_rap( true)
	gmg:set_smooth_on_surface_rim(false)

	-- gmg:set_damp(MinimalResiduumDamping())
	-- gmg:set_damp(0.8)
	-- gmg:set_damp(MinimalEnergyDamping())
	
	
	----------------------------------------------------------
	-- Linear Solver
	----------------------------------------------------------
	
		-- create Linear Solver
	GMresSolver = GMRES(20)
	GMresSolver:set_preconditioner(gmg)
	GMresSolver:set_convergence_check(LinearConvCheck)
	
	-- create Linear Solver
	BiCGStabSolver = BiCGStab()
	BiCGStabSolver:set_preconditioner(gmg)
	BiCGStabSolver:set_convergence_check(LinearConvCheck)

	gmgSolver = LinearSolver()
	gmgSolver:set_preconditioner(gmg)
	gmgSolver:set_convergence_check(LinearConvCheck)
	
	ilutSolver = LinearSolver()
	ilutSolver:set_preconditioner(ilu)
	ilutSolver:set_convergence_check(LinearConvCheck)
	
	
	-- choose a solver
	LinearSolver = BiCGStabSolver
	--LinearSolver = GMresSolver
	--LinearSolver = gmgSolver
	--LinearSolver = ilutSolver
	
	
	
	

	local NewtonSolverSteady = nil
	if self.doSteadyState then
		NewtonSolverSteady = NewtonSolver()
		NewtonSolverSteady:set_linear_solver(LinearSolver)
		NewtonSolverSteady:set_convergence_check(NewtonSteadyConvCheck)
		NewtonSolverSteady:set_line_search(NewtonLineSearch)
		
		if self.NewtonSteadyDebug then
			local dbgWriter_steady = GridFunctionDebugWriter(approxSpace)
			dbgWriter_steady:set_vtk_output(true)
			dbgWriter_steady:set_conn_viewer_output(false)
			dbgWriter_steady:set_base_dir(self.debug_dir)
			NewtonSolverSteady:set_debug(dbgWriter_steady)
		end
		
	end
		
	local limex = nil
	local NLSolver = NewtonSolver()
	
	
	if self.timeMethod == "limex" then
		NLSolver:set_linear_solver(LinearSolver)
		NLSolver:set_convergence_check(LimexConvCheck)
		NLSolver:auto_update(false)
		limex = myProblem:LimexObject( domainDisc, NLSolver)
	else
		NLSolver:set_linear_solver(LinearSolver)
		NLSolver:set_convergence_check(NewtonConvCheck)
		NLSolver:set_line_search(NewtonLineSearch)
		NLSolver:set_reassemble_J_freq(0)
		NLSolver:auto_update(true)
		op = AssembledOperator(timeDisc)
		op:init()
		NLSolver:init(op)
		if NLSolver:prepare(u) == false then
			print ("Newton solver prepare failed.") return op, NLSolver, NewtonSolverSteady, limex, 0
		end
	end
	
	if self.NewtonDebug then
		local dbgWriter = GridFunctionDebugWriter(approxSpace)
		dbgWriter:set_vtk_output(true)
		dbgWriter:set_conn_viewer_output(false)
		dbgWriter:set_base_dir(self.debug_dir)
		NLSolver:set_debug(dbgWriter)
					
	end
	
	TransientNewtonUpdater = nil
	if (self.NewtonUpdater and self.timeMethod == "euler")then
		TransientNewtonUpdater = NewtonUpdaterProjection()
		TransientNewtonUpdater:set_projection_fct(self.dim+1)
		TransientNewtonUpdater:set_max_threshold(1.1)
		TransientNewtonUpdater:set_min_threshold(-1e-05)
		NLSolver:setNewtonUpdater(TransientNewtonUpdater)
	end
	
	
	self.NewtonSolverDescSteady = NewtonSolverDescSteady
	self.NewtonSolverDesc = NewtonSolverDesc
	self.approxSpace = approxSpace
	self.boolSolverDesc = true
	
	  
	return op, NLSolver, NewtonSolverSteady, limex, 1
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
	return time_work_steady, linsolver_calls, linsolver_steps, 1.0

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
	local DTmax = self.DTmax
	local modifyDT = self.modifyDT
	local StepDebug = self.StepDebug
	local AbsDefect = self.AbsDefect
	local RedDefect = self.RedDefect
	local maxConvRate = self.maxConvRate
	local minConvRate = self.minConvRate

	if(self.boolSolverDesc == false) then print("SolverDesc Not initialized ") return 0, 0, 0, 0, 0, 0 end
	local solverDesc = self.NewtonSolverDesc
	local approxSpace = self.approxSpace
		
    local Newton_Steps = 0
    local Newton_Steps_fail = 0
    local linsolver_calls_step = 0
    local linsolver_steps_step = 0
    local CompletedStep = false
    local time = StartTime
    local time2 = StartTime
    local dt_in = dt
	if(modifyDT == false) then dt_in = DTmax end
        

        
	while CompletedStep==false  do
        Newton_Steps = Newton_Steps+1
		-- choose time step

		do_dt = math.min(dt_in,math.max((time+DTmax-time2), 0.0))
		print("Size of timestep dt: " .. do_dt)
		-- setup time Disc for old solutions and timestep
		timeDisc:prepare_step(solTimeSeries, do_dt)
	
		-- prepare newton solver
		if solver:prepare(u) == false then
			print ("Newton solver failed at step "..step.."."); return 0, 0, 0, 0, 0, 0;
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
				if StepDebug then
				
					convCheck =
						{
							type		= "standard",
							iterations	= 5,
							absolute	= self.AbsDefect,
							reduction	= self.RedDefect,
							verbose		= true
						}

					solver:set_convergence_check(ConvCheck(5, self.AbsDefect, self.RedDefect, true))
					local dbgWriter = GridFunctionDebugWriter(approxSpace)
					dbgWriter:set_vtk_output(true)
					dbgWriter:set_base_dir(self.debug_dir .."Steady")
					solver:set_debug(dbgWriter)
					solver:init(op)
					timeDisc:prepare_step(solTimeSeries, do_dt)
					if solver:prepare(u) == false then
						print ("Newton solver failed at DEBUG step "..step..".");
					end
					print("++++++ DEBUG STEP  BEGIN ++++++")
					solver:apply(u)
					print ("Newton solver failed at DEBUG step "..step..".");
				end
				
				return 0, 0, 0, 0, 0, 0;

			else
				VecScaleAssign(u, 1.0, solTimeSeries:latest())
			
			end
		else
            
            if  ((time2 + do_dt + (1e-07) * DTmin -time)/DTmax > 1.0 or modifyDT == false) then
            
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
            
            
            --[[if modifyDT then
                if(average_non_linear_rates>maxConvRate and false) then
                    dt_in = math.max(do_dt*red_factor_success,1.00001*DTmin)
                    print ("-------------------------------------------------------------------Time step decrease at Step " .. step-1 + (time2-time)/DTmax .. ", dt =  " .. dt_in .. ". ")
                else if (CompletedStep== false or Dt_factor < 0.98 ) then
						if(average_non_linear_rates<minConvRate or num_newton_steps<=optimal_newton_steps) then
							dt_in=math.min(incr_factor*dt_in,DTmax)
							print ("-------------------------------------------------------------------Time step increased at Step " .. step-1 + (time2-time)/DTmax ..", dt =  " .. dt_in .. ". ")
						end
                    end
                end
            end]]
            
			if modifyDT then
                if(num_newton_steps>30) then
                    dt_in = math.max(do_dt*red_factor_success,1.00001*DTmin)
                    print ("-------------------------------------------------------------------Time step decrease at Step " .. step-1 + (time2-time)/DTmax .. ", dt =  " .. dt_in .. ". ")
                else if (CompletedStep== false or Dt_factor < 0.98 ) then
						if(average_non_linear_rates<minConvRate or num_newton_steps<=optimal_newton_steps) then
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

  return Newton_Steps, Newton_Steps_fail, linsolver_calls_step, linsolver_steps_step, 1
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
	local dtlimex = self.DTmin

	--  Euclidean (algebraic) norm
	--local estimator = Norm2Estimator()



	--print (estimator)
	local limexEstimator = CompositeGridFunctionEstimator()
	limexEstimator:add(H1SemiComponentSpace("u", 2 ))
	limexEstimator:add(H1SemiComponentSpace("v", 2 ))

	--limexEstimator:add(L2ComponentSpace("u", 2))
	--limexEstimator:add(L2ComponentSpace("v", 2))
	

	limexEstimator:add(H1SemiComponentSpace("p", 2, ConstUserMatrix(0.0001) ))
	--limexEstimator:add(L2ComponentSpace("p", 2))
	limexEstimator:add(L2ComponentSpace("c", 2))
	
	--limexEstimator:use_strict_relative_norms(1)
	print(limexEstimator:config_string())

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
	  rhoSafetyOPT = 0.25,
	  dtred = 0.5,
	  dtIncr = 1.5,
	  matrixCache = false,
	  conservative = true
	  
	}


	-- setup for time integrator
	local limex = util.limex.CreateIntegrator(limexDesc)

	limex:set_time_step(limexDesc.dt)
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
	
	limex:set_start_step(1)
	limex:apply(u, EndTime, u, StartTime)
	n_step = limex:get_step()-1
	
    local Newton_Steps = NLSolver:total_linsolver_calls()
	local Newton_Steps_fail = 0
	linsolver_calls_step = NLSolver:total_linsolver_calls()
	linsolver_steps_step =  NLSolver:total_linsolver_steps()
	limex:set_time_step(self.DTmax/n_step)
	NLSolver:clear_average_convergence();

  return Newton_Steps, Newton_Steps_fail, linsolver_calls_step, linsolver_steps_step, 1
end

return myProblem
    
