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
	self.NewtonDebug = o.NewtonDebug
	self.NewtonSteadyDebug = o.NewtonSteadyDebug
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
	self.AbsDefect_steady = 1e-08
	self.RedDefect_steady = 1e-05

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
	self.boolLoadCheckPoint = o.boolLoadCheckPoint
	self.boolSaveCheckPoint = o.boolSaveCheckPoint
	
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
	--local NewtonSteadyConvCheck=ConvCheck(self.max_newton_steps_steady_state, self.AbsDefect_steady, self.RedDefect_steady, true)
	--local NewtonConvCheck=ConvCheck(self.max_newton_steps_transient, self.AbsDefect, self.RedDefect, true)
	
	local NewtonSteadyConvCheck = CompositeConvCheck(approxSpace)
	NewtonSteadyConvCheck:set_maximum_steps(50)
	NewtonSteadyConvCheck:set_group_check({"u", "v"}, 1.0e-12, 1.0e-8)
	NewtonSteadyConvCheck:set_component_check("p", 1.0e-12, 1.0e-8)
	NewtonSteadyConvCheck:set_component_check("c", 1.0e-12, 1.0e-8)
	NewtonSteadyConvCheck:set_component_check("k", 1.0e-12, 1.0e-8)
	NewtonSteadyConvCheck:set_component_check("omega", 1.0e-12, 1.0e-8)
	NewtonSteadyConvCheck:disable_rest_check()
	NewtonSteadyConvCheck:set_verbose(true)
	
	local NewtonConvCheck = CompositeConvCheck(approxSpace)
	NewtonConvCheck:set_maximum_steps(self.max_newton_steps_transient)
	NewtonConvCheck:set_group_check({"u", "v"}, self.AbsDefect, self.RedDefect)
	NewtonConvCheck:set_component_check("p", self.AbsDefect, self.RedDefect)
	NewtonConvCheck:set_component_check("c", self.AbsDefect, self.RedDefect)
	NewtonConvCheck:set_component_check("k", self.AbsDefect, self.RedDefect)
	NewtonConvCheck:set_component_check("omega", self.AbsDefect, self.RedDefect)
	NewtonConvCheck:disable_rest_check()
	NewtonConvCheck:set_verbose(true)
	
	
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
	--ilu:set_ordering_algorithm(TopologicalOrdering())
	--ilu:set_sort(true)
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
	gmg:set_num_presmooth(0)
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
	--LinearSolver=baseSolver
	
	self.LinearSolver = LinearSolver
	
	

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



	TransientNewtonUpdater = NewtonUpdaterProjection()
	TransientNewtonUpdater:set_projection_fct(self.dim+1)
	TransientNewtonUpdater:set_max_threshold(2.0)
	TransientNewtonUpdater:set_min_threshold(-2.0)
	
		
	local limex = nil
	local NLSolver = NewtonSolver()
	
	
	if self.timeMethod == "limex" then
		NLSolver:set_linear_solver(LinearSolver)
		NLSolver:set_convergence_check(LimexConvCheck)
		--NLSolver:setNewtonUpdater(TransientNewtonUpdater)
		limex = myProblem:LimexObject( domainDisc, NLSolver)
	else
		NLSolver:set_linear_solver(LinearSolver)
		NLSolver:set_convergence_check(NewtonConvCheck)
		NLSolver:set_line_search(NewtonLineSearch)
		NLSolver:set_reassemble_J_freq(0)
		--NLSolver:setNewtonUpdater(TransientNewtonUpdater)
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
	local fixer = DirichletBoundary()
	domainDisc:add(fixer)
	fixer:invert_subset_selection()
	fixer:add("c", "")
	fixer:add("k", "")
	fixer:add("omega", "")
	

	solver:init(AssembledOperator(domainDisc))
	
	solver:prepare(u)
	
	domainDisc:adjust_solution(u)
	KinTurbulentViscosity:update()

	-- apply the solver for the stationary pressure problem
     print("++++++ STEADY STATE CALCULATION BEGIN ++++++")
	tBefore_s= os.clock()
	if not solver:apply(u) then
		print("===> THE PREPARATION PHASE FAILED! <===")
		domainDisc:remove (fixer)
		return 0, 0, 0, 0
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
	local DTmax = EndTime-StartTime
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

  return Newton_Steps, Newton_Steps_fail, linsolver_calls_step, linsolver_steps_step, dt_in, 1
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

	limexEstimator:add(H1SemiComponentSpace("p", 2))--, ConstUserMatrix(1/207414416) ))
	limexEstimator:add(L2ComponentSpace("c", 2))
	
	limexEstimator:use_strict_relative_norms(1)

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
	  matrixCache = true,
	  conservative = false
	  
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

	limex:apply(u, EndTime, u, StartTime)
	
    local Newton_Steps = NLSolver:total_linsolver_calls()
	local Newton_Steps_fail = 0
	linsolver_calls_step = NLSolver:total_linsolver_calls()
	linsolver_steps_step =  NLSolver:total_linsolver_steps()
	NLSolver:clear_average_convergence();

  return Newton_Steps, Newton_Steps_fail, linsolver_calls_step, linsolver_steps_step, 1
end

myProblem.CheckJacobian = function(self, domainDisc, u, approxSpace, perturbFct, solveLinearSystem)

	------------------------------------------------------------
	-- Settings
	------------------------------------------------------------

	if solveLinearSystem == nil then
		solveLinearSystem = true
	end

	local linSolver = self.LinearSolver

	local innerSubsets = table.concat(Inner_total, ",")
	local boundarySubsets = "Inlet,UpperWall,LowerWall,CylinderWall,Outlet"

	local epsList = {
		1.0e-3,
		1.0e-4,
		1.0e-5,
		1.0e-6,
		1.0e-7,
		1.0e-8
	}


	------------------------------------------------------------
	-- UPDATE STATE-DEPENDENT PARAMETERS
	--
	-- Add here all quantities that must be updated whenever
	-- the current solution changes.
	------------------------------------------------------------

	local function UpdateParameters()
		self.KinTurbulentViscosity:update()
	end


	------------------------------------------------------------
	-- Jacobian perturbation function
	--
	-- Must be GLOBAL because UG4 Interpolate searches for the
	-- callback by its string name in the global Lua namespace.
	------------------------------------------------------------

	function JacobianPerturbation(x, y, t)
		return x
	end


	------------------------------------------------------------
	-- Build perturbation direction v
	------------------------------------------------------------

	local v = u:clone()
	v:set(0.0)

	Interpolate("JacobianPerturbation", v, perturbFct, innerSubsets)
	Interpolate(0.0, v, perturbFct, boundarySubsets)


	------------------------------------------------------------
	-- Working vectors
	------------------------------------------------------------

	local rNewton = u:clone()
	rNewton:set(0.0)

	local r0 = u:clone()
	r0:set(0.0)

	local r1 = u:clone()
	r1:set(0.0)

	local Jv = u:clone()
	Jv:set(0.0)

	local fd = u:clone()
	fd:set(0.0)

	local err = u:clone()
	err:set(0.0)


	------------------------------------------------------------
	-- Adjust current solution
	------------------------------------------------------------

	domainDisc:adjust_solution(u)


	------------------------------------------------------------
	-- Newton defect
	--
	-- Newton computes the defect BEFORE the step update.
	------------------------------------------------------------

	domainDisc:assemble_defect(rNewton, u)


	------------------------------------------------------------
	-- Update state-dependent parameters
	--
	-- This corresponds to the Newton step update before the
	-- Jacobian is assembled.
	------------------------------------------------------------

	UpdateParameters()


	------------------------------------------------------------
	-- Assemble Jacobian
	------------------------------------------------------------

	local J = AssembledLinearOperator(domainDisc)
	domainDisc:assemble_jacobian(J, u)


	------------------------------------------------------------
	-- Assemble the reference defect for the FD test
	--
	-- IMPORTANT:
	-- This is assembled AFTER UpdateParameters(), so that
	--
	--     R(u), J(u), R(u + eps*v)
	--
	-- all use exactly the same frozen parameters.
	------------------------------------------------------------

	r0:set(0.0)
	domainDisc:assemble_defect(r0, u)


	------------------------------------------------------------
	-- Component spaces
	------------------------------------------------------------

	local spaceUV = GridFunctionComponentSpace("u,v")
	local spaceP = GridFunctionComponentSpace("p")
	local spaceC = GridFunctionComponentSpace("c")
	local spaceK = GridFunctionComponentSpace("k")
	local spaceOmega = GridFunctionComponentSpace("omega")

	local spaceKInner = GridFunctionComponentSpace("k", innerSubsets)
	local spaceOmegaInner = GridFunctionComponentSpace("omega", innerSubsets)

	local spacePertInner = GridFunctionComponentSpace(perturbFct, innerSubsets)
	local spacePertInlet = GridFunctionComponentSpace(perturbFct, "Inlet")
	local spacePertOutlet = GridFunctionComponentSpace(perturbFct, "Outlet")
	local spacePertUpper = GridFunctionComponentSpace(perturbFct, "UpperWall")
	local spacePertLower = GridFunctionComponentSpace(perturbFct, "LowerWall")
	local spacePertCylinder = GridFunctionComponentSpace(perturbFct, "CylinderWall")


	------------------------------------------------------------
	-- Optional Newton / linear-system analysis
	------------------------------------------------------------

	if solveLinearSystem then

		--------------------------------------------------------
		-- Actual Newton correction
		--
		-- J * corr = R
		-- u_new = u - lambda * corr
		--------------------------------------------------------

		local corr = u:clone()
		corr:set(0.0)

		linSolver:init(J, corr)

		if not linSolver:apply(corr, rNewton) then

			print("")
			print("================================================")
			print(" LINEAR SOLVE FOR NEWTON CORRECTION FAILED")
			print("================================================")
			print("")

		else

			print("")
			print("================================================")
			print(" ACTUAL NEWTON CORRECTION")
			print("================================================")

			print(string.format("|corr|         = %.6e", VecNorm(corr)))
			print(string.format("|corr u,v|     = %.6e", spaceUV:norm(corr)))
			print(string.format("|corr p|       = %.6e", spaceP:norm(corr)))
			print(string.format("|corr c|       = %.6e", spaceC:norm(corr)))
			print(string.format("|corr k|       = %.6e", spaceK:norm(corr)))
			print(string.format("|corr omega|   = %.6e", spaceOmega:norm(corr)))
			print(string.format("|corr k inner| = %.6e", spaceKInner:norm(corr)))
			print(string.format("|corr w inner| = %.6e", spaceOmegaInner:norm(corr)))

			print("================================================")
			print("")


			----------------------------------------------------
			-- Norm calls may change parallel storage.
			----------------------------------------------------

			corr:enforce_consistent_type()
			u:enforce_consistent_type()


			----------------------------------------------------
			-- Save original solution
			----------------------------------------------------

			local uSave = u:clone()

			corr:enforce_consistent_type()
			uSave:enforce_consistent_type()


			----------------------------------------------------
			-- Newton line-search diagnostic
			----------------------------------------------------

			local lambdaList = {
				1.0,
				0.5,
				0.25,
				0.125,
				0.0625,
				0.03125,
				0.015625,
				0.0078125,
				0.00390625,
				0.001953125,
				0.0009765625,
				0.00048828125,
				0.000244140625,
				0.0001220703125
			}

			print("NEWTON TRIAL DEFECTS")

			for _, lambda in ipairs(lambdaList) do

				corr:enforce_consistent_type()
				uSave:enforce_consistent_type()
				u:enforce_consistent_type()

				VecScaleAdd2(u, 1.0, uSave, -lambda, corr)

				domainDisc:adjust_solution(u)

				------------------------------------------------
				-- DO NOT call UpdateParameters() here.
				--
				-- The Newton line search uses the parameters
				-- frozen at the step update.
				------------------------------------------------

				local rTrial = u:clone()
				rTrial:set(0.0)

				domainDisc:assemble_defect(rTrial, u)

				print(string.format(
					"lambda = %.6f   |R(u-lambda*corr)| = %.6e",
					lambda,
					VecNorm(rTrial)
				))

				u:enforce_consistent_type()
				uSave:enforce_consistent_type()

				VecAssign(u, uSave)
			end


			----------------------------------------------------
			-- Restore original solution and parameter state
			----------------------------------------------------

			u:enforce_consistent_type()
			uSave:enforce_consistent_type()

			VecAssign(u, uSave)

			domainDisc:adjust_solution(u)

			UpdateParameters()

			print("")
		end

	else

		print("")
		print("================================================")
		print(" LINEAR SYSTEM / NEWTON ANALYSIS DISABLED")
		print("================================================")
		print("")

	end


	------------------------------------------------------------
	-- Re-establish state for Jacobian finite-difference test
	------------------------------------------------------------

	domainDisc:adjust_solution(u)
	UpdateParameters()

	------------------------------------------------------------
	-- Reassemble reference defect with the same frozen
	-- parameters used for the FD perturbations.
	------------------------------------------------------------

	r0:set(0.0)
	domainDisc:assemble_defect(r0, u)


	------------------------------------------------------------
	-- Compute J*v
	------------------------------------------------------------

	Jv:set(0.0)

	v:enforce_consistent_type()
	Jv:enforce_consistent_type()

	J:apply(Jv, v)


	------------------------------------------------------------
	-- Print perturbation information
	------------------------------------------------------------

	print("")
	print("================================================")
	print(" " .. string.upper(perturbFct) .. "-DIRECTION JACOBIAN FINITE-DIFFERENCE TEST")
	print("================================================")
	print("")

	print(string.upper(perturbFct) .. " PERTURBATION v ON SUBSETS")

	print(string.format("inner    : %.6e", spacePertInner:norm(v)))
	print(string.format("inlet    : %.6e", spacePertInlet:norm(v)))
	print(string.format("outlet   : %.6e", spacePertOutlet:norm(v)))
	print(string.format("upper    : %.6e", spacePertUpper:norm(v)))
	print(string.format("lower    : %.6e", spacePertLower:norm(v)))
	print(string.format("cylinder : %.6e", spacePertCylinder:norm(v)))

	print("")


	------------------------------------------------------------
	-- Restore consistent storage after norm calls
	------------------------------------------------------------

	v:enforce_consistent_type()
	u:enforce_consistent_type()
	r0:enforce_consistent_type()
	Jv:enforce_consistent_type()


	------------------------------------------------------------
	-- Finite-difference Jacobian test
	------------------------------------------------------------

	for _, eps in ipairs(epsList) do

		local uPert = u:clone()

		u:enforce_consistent_type()
		v:enforce_consistent_type()
		uPert:enforce_consistent_type()

		VecScaleAdd2(uPert, 1.0, u, eps, v)

		domainDisc:adjust_solution(uPert)


		--------------------------------------------------------
		-- IMPORTANT:
		--
		-- Do NOT call UpdateParameters() here.
		--
		-- The state-dependent parameters remain frozen at the
		-- state at which J was assembled.
		--------------------------------------------------------

		r1:set(0.0)
		domainDisc:assemble_defect(r1, uPert)


		--------------------------------------------------------
		-- FD = [R(u + eps*v) - R(u)] / eps
		--------------------------------------------------------

		r1:enforce_consistent_type()
		r0:enforce_consistent_type()
		fd:enforce_consistent_type()

		VecScaleAdd2(fd, 1.0 / eps, r1, -1.0 / eps, r0)


		--------------------------------------------------------
		-- err = FD - J*v
		--------------------------------------------------------

		fd:enforce_consistent_type()
		Jv:enforce_consistent_type()
		err:enforce_consistent_type()

		VecScaleAdd2(err, 1.0, fd, -1.0, Jv)


		--------------------------------------------------------
		-- Global norms
		--------------------------------------------------------

		local normFD = VecNorm(fd)
		local normJv = VecNorm(Jv)
		local normErr = VecNorm(err)

		local denom = math.max(normFD, normJv, 1.0e-30)
		local relErr = normErr / denom


		--------------------------------------------------------
		-- Detailed block analysis at eps = 1e-6
		--------------------------------------------------------

		if eps == 1.0e-6 then

			print("")
			print("COMPONENT-WISE JACOBIAN CHECK AT eps = 1e-6")
			print("Input block: " .. perturbFct)
			print("")

			print(string.format(
				"u,v   : |FD| = %.6e   |Jv| = %.6e   |err| = %.6e",
				spaceUV:norm(fd),
				spaceUV:norm(Jv),
				spaceUV:norm(err)
			))

			print(string.format(
				"p     : |FD| = %.6e   |Jv| = %.6e   |err| = %.6e",
				spaceP:norm(fd),
				spaceP:norm(Jv),
				spaceP:norm(err)
			))

			print(string.format(
				"c     : |FD| = %.6e   |Jv| = %.6e   |err| = %.6e",
				spaceC:norm(fd),
				spaceC:norm(Jv),
				spaceC:norm(err)
			))

			print(string.format(
				"k     : |FD| = %.6e   |Jv| = %.6e   |err| = %.6e",
				spaceK:norm(fd),
				spaceK:norm(Jv),
				spaceK:norm(err)
			))

			print(string.format(
				"omega : |FD| = %.6e   |Jv| = %.6e   |err| = %.6e",
				spaceOmega:norm(fd),
				spaceOmega:norm(Jv),
				spaceOmega:norm(err)
			))

			print("")

			print(string.format(
				"k inner     : |FD| = %.6e   |Jv| = %.6e   |err| = %.6e",
				spaceKInner:norm(fd),
				spaceKInner:norm(Jv),
				spaceKInner:norm(err)
			))

			print(string.format(
				"omega inner : |FD| = %.6e   |Jv| = %.6e   |err| = %.6e",
				spaceOmegaInner:norm(fd),
				spaceOmegaInner:norm(Jv),
				spaceOmegaInner:norm(err)
			))

			print("")
		end


		--------------------------------------------------------
		-- Global FD comparison
		--------------------------------------------------------

		print(string.format(
			"eps = %.1e   |FD| = %.6e   |Jv| = %.6e   |FD-Jv| = %.6e   rel = %.6e",
			eps,
			normFD,
			normJv,
			normErr,
			relErr
		))


		--------------------------------------------------------
		-- Restore compatible storage for next epsilon
		--------------------------------------------------------

		u:enforce_consistent_type()
		v:enforce_consistent_type()
		r0:enforce_consistent_type()
		r1:enforce_consistent_type()
		Jv:enforce_consistent_type()
		fd:enforce_consistent_type()
		err:enforce_consistent_type()
	end


	print("================================================")
	print("")

end

--------------------------------------------------------------------------------
-- CheckPoint
--------------------------------------------------------------------------------
myProblem.SaveCheckPoint = function (self,u,folder)

	if (not(self.boolSaveCheckPoint)) then return end
	
	SaveToFile(u, folder .. "/CheckPoint" .. ".vec")
	print("Saving CheckPoint: DONE")

end

myProblem.LoadCheckPoint = function (self,u,folder)

	
	local boolInterpolate = false
	local filename = folder .. "/Integral.txt"
	local ChekPointname = folder .. "/CheckPoint.vec"
	
	local step = 0
	local time_work_total = 0.0
	local Checkfile = io.open(ChekPointname, "r")

	if Checkfile then
		print("Loading CheckPoint")
		Checkfile:close()
		ReadFromFile(u, ChekPointname)
		

	else
		print("File does not exist: " .. ChekPointname)
		boolInterpolate = true
	end

	
	return time, step, time_work_total, boolInterpolate
end

return myProblem
    
