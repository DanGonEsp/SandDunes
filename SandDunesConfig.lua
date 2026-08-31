--------------------------------------------------------------------------------
-- Initial Parameters
--------------------------------------------------------------------------------
local myProblem = {}

myProblem.Init = function(self, o)

		-- Numerical parameters of the discretization
	self.dim = o.dim
	self.dir_name = o.dir_name
	self.file_name = o.file_name
	self.folder_name = o.folder_name
	self.elem_type = o.elem_type
	self.numRefs = o.numRefs
	self.numPreRefs = o.numPreRefs
	
	self.simCase = o.simCase
	
	self.timeMethod  = o.timeMethod
	self.DT = o.DT
	self.DTmax = o.DTmax
	self.DTmin = o.DTmin
	self.DTLimex = o.DTLimex
	self.numTimeSteps = o.numTimeSteps
	self.startTime = o.startTime
	self.endTime = o.endTime
	self.modifyDT  = o.modifyDT
	
	
		--Output Data
	self.boolData = o.boolData
	self.data_name = o.data_name
	self.outputFactor = o.outputFactor
	self.boolCheckPoint = o.boolCheckPoint

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
	self.VelErrorNorm = o.VelErrorNorm
	self.PressErrorNorm = o.PressErrorNorm
	self.VolErrorNorm = o.VolErrorNorm
	self.alphaVel  = o.alphaVel
	self.alphaPress = o.alphaPress
	self.alphaVol = o.alphaVol

	self.max_newton_steps_steady_state = o.max_newton_steps_steady_state
	self.max_newton_steps_transient = o.max_newton_steps_transient
	
	self.SteadyAbsDefect = o.SteadyAbsDefect
	self.SteadyRedDefect = o.SteadyRedDefect
	
	self.AbsDefect = o.AbsDefect
	self.RedDefect = o.RedDefect

	self.LinAbsDefectImp = o.LinAbsDefectImp
	self.LinRedDefectImp = o.LinRedDefectImp
	self.LinAbsDefectLim = o.LinAbsDefectLim
	self.LinRedDefectLim = o.LinRedDefectLim
	self.max_linear_steps_Lim = o.max_linear_steps_Lim
	self.max_linear_steps_Imp = o.max_linear_steps_Imp
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
	self.boolSlipVel = o.boolSlipVel
	self.boolpress_jump = o.boolpress_jump
	self.boolNormal = o.boolNormal
	self.boolFixVel = o.boolFixVel
	self.boolFixVol = o.boolFixVol
	self.boolMassTerm = o.boolMassTerm
	self.boolDensityMean = o.boolDensityMean
	
	self.inflow = o.inflow
	self.SlipVelValue = o.SlipVelValue
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
	self.div_correction = o.div_correction
	self.boolIPVelocity = o.boolIPVelocity
	self.boolTransportJac = o.boolTransportJac
	self.turbViscMethod = o.turbViscMethod
	self.modellconstant = o.modellconstant
	self.update_turb = o.update_turb

	
	--Material Properties
	self.nu_a = o.nu_a
	self.rho_a = o.rho_a
	self.rho_s = o.rho_s
	self.dp = o.dp
	self.nu_s = o.nu_s
	self.c_init = o.c_init

	self.alpha_max = o.alpha_max
	self.alpha_min = o.alpha_min
	self.packing_factor = o.packing_factor
	self.grad_limit = o.grad_limit
	self.slope_limit = o.slope_limit
	self.granular_model = o.granular_model
	self.density_model = o.density_model
	self.interface_value = o.interface_value
	self.drag_mod = o.drag_mod
	self.riemman = o.riemman


	self.FR = o.FR
	self.B_phi = o.B_phi
	self.deltaGamma = o.deltaGamma
	self.Visc_limit = o.Visc_limit

	self.deltaPs = o.deltaPs
	self.deltaI = o.deltaI
	self.FricMu_1 = o.FricMu_1
	self.FricMu_2 = o.FricMu_2
	self.I_0 = o.I_0
	self.gravity = o.gravity
	
	
	self.boolSolverDesc = false
	self.NewtonSolverDescSteady = {}
	self.NewtonSolverDesc = {}
	
	
	self.gridName = o.gridName
  
  
end
--------------------------------------------------------------------------------
-- File Names
--------------------------------------------------------------------------------
myProblem.FileNames = function (self,rank,SpaceSize)


	riemman_name = nil
	if self.riemman == 0 then
		riemman_name = "Upwind"
	elseif self.riemman == 1 then
		riemman_name = "Godunov"
	elseif self.riemman == 2 then
		riemman_name = "Rusanov"
	elseif self.riemman == 3 then
		riemman_name = "Roe"
	else
		print ("Numerical Flux Scheme for NonLinear Scalar conservation Law not defined"); exit();
	end

	local folder_name = nil
	
	if self.dir_name == "" then
		folder_name = self.folder_name
	else
		folder_name = self.dir_name .. "/" .. self.folder_name
	end

	---------------------
	-- FolderName
	---------------------

	folder_name = folder_name .. "_" .. self.dim.. "D"
	folder_name = folder_name .. "-" .. self.timeMethod
	
	if self.bStokes then
		folder_name = folder_name .. "-Stokes"
	end

	folder_name = folder_name .. "-" ..riemman_name

	if self.boolRelativeVel then
		folder_name =folder_name .. "-RelVel"
	else
		folder_name =folder_name .. "-NoRelVel"
	end


	if  self.boolGradientPsSource  or self.boolSource  then

		if self.boolSource then
			if self.consistentRho_in_source then
				folder_name =folder_name .. "-Consistent"
			end
			folder_name = folder_name .. "-MG_Force"
		end
		
		if self.boolGradientPsSource then
			folder_name = folder_name .. "-DPs"
		else
			folder_name = folder_name .. "-DPh"
		end
		
	else
		folder_name = folder_name .. "-NoForce"
	end

	if self.boolAveDiff then
		folder_name = folder_name .. "-AveDiff"
	else
		folder_name = folder_name .. "-NoAveDiff"
	end

	if (self.boolSlipDiff and self.boolSlipVel) then print ("SlipVel and SlipDiff activated at the same time."); exit(); end
	if self.boolSlipDiff then
		folder_name = folder_name .. "-SlipDiff"
	else if self.boolSlipVel then
			folder_name = folder_name .. "-SlipVel"
		else
			folder_name = folder_name .. "-NoSlip"
		end
	end
	

	if (rank==0) then
		if not DirectoryExists (folder_name) then
			CreateDirectory (folder_name)
		end
		print(folder_name)
	end
	
	---------------------
	-- VTK folder
	---------------------
	local file_name = self.file_name
	
	vtk_file_name = file_name .. "-" .. self.elem_type
	vtk_file_name = vtk_file_name .. "-lev" .. self.numRefs
	
	local fixedSpaceSize = string.format("%04d", SpaceSize)
	vtk_file_name = vtk_file_name .. "_Parallel_" .. fixedSpaceSize

	folder_vtk = folder_name .. "/" .. vtk_file_name
	
	if (rank== 0) then
		
		if not DirectoryExists (folder_vtk) then
			CreateDirectory (folder_vtk)
		else
			if(self.NewtonDebug or self.NewtonSteadyDebug) then
				local cmd = "rm -rf " .. folder_vtk
				os.execute(cmd)
				print(" Directory reseted")
				CreateDirectory (folder_vtk)
			end
		end
	end

	

	vtk_file_name = folder_vtk .. "/" .. self.file_name -- VTK output file name base
	
	
	self.debug_dir = folder_vtk
	self.riemman_name = riemman_name

	print("Files Setting DONE")
	print("File Adress: ".. vtk_file_name)
	
	-- vtk_file_name: Name of the vtk files
	-- folder_vtk: Directory where vtk files are saved
	-- folder_name: Directory where the different simulations are saved, e.g. folder2dLimex, folder3dEuler
	
	return vtk_file_name,folder_vtk,folder_name
end

--------------------------------------------------------------------------------
-- LOG Files
--------------------------------------------------------------------------------
myProblem.LogFiles = function (self,rank_t,Name)
	
	
	local counter = 0
	local filename = Name .. "_" .. counter
	while io.open(filename, "r") do

		counter = counter + 1
		filename = Name .. "_" .. counter

	end
	GetLogAssistant():enable_file_output(true, filename)
	if rank_t > 0 then GetLogAssistant():enable_terminal_output(false) end

end





--------------------------------------------------------------------------------
-- Printing Simulation Settings
--------------------------------------------------------------------------------
myProblem.PrintingSettings = function (self)
	-----------------------------------------------------------------------
	print (" Sand Dune Dynamics     " .. os.date("%A, %B %d, %Y at %I:%M %p"))
	print (" Geometry: " .. self.gridName ..", dim = " .. self.dim)
	print (" Physical parameter:")
	print ("	Table case		= " .. self.simCase + 1)
	print ("	inflow			= " .. self.inflow)
	print ("	SlipVel			= " .. self.SlipVelValue)
	print ("	Stokes			= " .. tostring (self.bStokes))
	print ("	Steady state    	= " .. tostring (self.doSteadyState))
	print ("	BodyForce       	= " .. tostring (self.boolSource))
	print ("	Consisten Gravity	= " .. tostring (self.consistentRho_in_source))
	print ("	Relative Vel    	= " .. tostring (self.boolRelativeVel))
	print ("	Ps gradient     	= " .. tostring (self.boolGradientPsSource))
	print ("	Ps in visc      	= " .. tostring (self.boolViscPs))
	print ("	Diffusion       	= " .. tostring (self.boolAveDiff))
	print ("	SlipDiff         	= " .. tostring (self.boolSlipDiff))
	print ("	SlipVel         	= " .. tostring (self.boolSlipVel))
	print ("	MassMean         	= " .. tostring (self.boolDensityMean))
	-----------------------------------------------------------------------
	print (" Numerical parameter:")
	print ("	elem_type		= " .. self.elem_type)
	print ("	numRefs			= " .. self.numRefs)
	print ("	numPreRefs		= " .. self.numPreRefs)
	print ("	timeMethod		= " .. self.timeMethod)
	print ("	DT   			= " .. self.DT)
	print ("	numTimeSteps		= " .. self.numTimeSteps)
	print ("	noLaplace		= " .. tostring (self.bNoLaplace))
	print ("	exactJac		= " .. tostring (self.bExactJac))
	print ("	PecletBlend		= " .. tostring (self.bPecletBlend))
	print ("	upwind_m		= " .. self.upwind_m)
	print ("	upwind_t		= " .. self.upwind_t)
	print ("	Num Flux Scheme		= " .. self.riemman_name)
	print ("	pac			= " .. tostring (self.bPac))
	print ("	stab			= " .. self.stab)
	print ("	difflength		= " .. self.diffLength)
	print ("	Turbulence		= " .. self.turbViscMethod)
	-----------------------------------------------------------------------
	print (" Linear Solver parameters:")
	print ("	damping_mg		= " .. self.damping_mg)
	print ("	beta			= " .. self.value_beta)
	print ("	AbsDefImp		= " .. self.LinAbsDefectImp)
	print ("	RedDefImp		= " .. self.LinRedDefectImp)
	print ("	AbsDefLim 		= " .. self.LinAbsDefectLim)
	print ("	RedDefLim		= " .. self.LinRedDefectLim)
	print ("	MaxStepsLim		= " .. tostring (self.max_linear_steps_Lim))
	print ("	MaxStepsImp		= " .. tostring (self.max_linear_steps_Imp))

	-----------------------------------------------------------------------

	if self.timeMethod == "limex" then
		print (" Limex Numerical parameter:")
		print ("	Vel norm		= " .. self.VelErrorNorm)
		print ("	Press norm		= " .. self.PressErrorNorm)
		print ("	Vol  norm		= " .. self.VolErrorNorm)
		print ("	Vel norm scale		= " .. self.alphaVel)
		print ("	Press norm scale	= " .. self.alphaPress)
		print ("	Vol norm scale		= " .. self.alphaVol)
	end

end


--------------------------------------------------------------------------------
-- Approximation Space
--------------------------------------------------------------------------------
myProblem.ApproximationSpace = function (self,allSubsets)
	print("Approximation Space Setting")
	local fct_cmp_tbl = nil
	local vel_cmp_tbl = nil
			
	if self.dim == 3 then
		fct_cmp_tbl = {"u", "v", "w", "p", "c"}
		vel_cmp_tbl = {"u", "v", "w"}
	else
		fct_cmp_tbl = {"u", "v", "p", "c"}
		vel_cmp_tbl = {"u", "v"}
	end

	-- Create the domain, load the grid and refine it
	local dom = util.CreateDomain (self.gridName, self.numPreRefs)
	balancer.RefineAndRebalanceDomain (dom, self.numRefs - self.numPreRefs)
	--util.refinement.CreateRegularHierarchy(dom, self.numRefs, true)

	print ("Domain-info:")
	print (dom:domain_info():to_string())

	-- Create the vertex-centered approximation space
	local approxSpace = ApproximationSpace (dom)

	approxSpace:add_fct("u", "Lagrange",1,allSubsets)
	approxSpace:add_fct("v", "Lagrange",1,allSubsets)
	if self.dim == 3 then
		approxSpace:add_fct("w", "Lagrange",1,allSubsets)
	end
	approxSpace:add_fct("p", "Lagrange",1,allSubsets)
	approxSpace:add_fct("c", "Lagrange",1,allSubsets)

	approxSpace:init_levels()
	approxSpace:init_top_surface()
	approxSpace:print_statistic()

	OrderLex (approxSpace, "y")
	--OrderCuthillMcKee(approxSpace,true)

	util.solver.defaults.approxSpace = approxSpace


	-- grid function for the solution
	local u = GridFunction (approxSpace)
	u:set(0)

	self.fct_cmp_tbl = fct_cmp_tbl
	self.vel_cmp_tbl = vel_cmp_tbl
	self.u = u
	self.approxSpace = approxSpace
	
	print("Approximation Space DONE")
	
	
	return approxSpace,u
end

--------------------------------------------------------------------------------
-- Interface Parameters List
--------------------------------------------------------------------------------
myProblem.InterfaceParameters = function (self)
	
	local InterfaceValues = Interface()

	InterfaceValues:set_particle_diameter(self.dp)
	InterfaceValues:set_particle_density(self.rho_s)
	InterfaceValues:set_air_density(self.rho_a)
	InterfaceValues:set_fluid_Visc(self.nu_a*self.rho_a)
	InterfaceValues:set_particle_kinVisc(self.nu_s)
	InterfaceValues:set_alpha_max(self.alpha_max)
	InterfaceValues:set_alpha_min(self.alpha_min)
	InterfaceValues:set_packing_factor(self.packing_factor)
	InterfaceValues:set_FrictionMu_1(self.FricMu_1)
	InterfaceValues:set_FrictionMu_2(self.FricMu_2)
	InterfaceValues:set_I_0(self.I_0)
	InterfaceValues:set_deltaI(self.deltaI)
	InterfaceValues:set_deltaPs(self.deltaPs)
	InterfaceValues:set_FR(self.FR)
	InterfaceValues:set_B_phi(self.B_phi)
	InterfaceValues:set_deltaGamma(self.deltaGamma)
	InterfaceValues:set_limit(self.Visc_limit)
	InterfaceValues:set_bool_particle_pressure_force(self.boolGradientPsSource)
	InterfaceValues:set_bool_consistent_gravity(false)
	InterfaceValues:set_reference_pressure(self.ReferencePressure)
	InterfaceValues:set_time_step_factor(self.DT)
	InterfaceValues:set_interface_volume_fraction(self.interface_value)
	InterfaceValues:set_drag_model(self.drag_mod)
	InterfaceValues:set_relative_vel_error(1)
	InterfaceValues:set_bool_initialized(true)
	
	self.InterfaceValues = InterfaceValues
	
	print("Interfce Parameters list: DONE")
	
	return InterfaceValues
end

--------------------------------------------------------------------------------
-- Variables
--------------------------------------------------------------------------------
myProblem.Clousures = function (self,approxSpace,u,walls)
		
	local InterfaceValues = self.InterfaceValues
	-------------------------------------------------------------- VelocityGradMag

	local gamma = ShearStressFV1(approxSpace,u)
	self.gamma = gamma

	---------------------------------------------------------------------- Density
	local Density = nil
	if self.bStokes then
		Density = ConstUserNumber(self.rho_a)
	else
		Density = GranularDensityLinker();
		Density:set_model(self.density_model)
		Density:set_phase_parameters(self.InterfaceValues)
	end


	self.Density = Density

	---------------------------------------------------------------------- Viscosity

	local Inverse_RHO = InverseLinker();
	Inverse_RHO:divide(1.0,Density);

	Scale_RHO = ScaleLinker();
	Scale_RHO:set_import_1(Inverse_RHO)
	Scale_RHO:set_import_2(self.rho_a)

	local KinMixViscosity = nil
	KinMixViscosity = ScaleLinker();
	KinMixViscosity:set_import_1(Inverse_RHO)
	

	if self.turbViscMethod=="dyn" then
		KinTurbulentViscosity = FV1DynamicTurbViscData(approxSpace,u)
		
	else
		KinTurbulentViscosity = FV1SmagorinskyTurbViscData(approxSpace,u,self.modellconstant)
	 end
	KinTurbulentViscosity:set_turbulence_zero_bnd(walls)
	KinTurbulentViscosity:set_kinematic_viscosity(0.0)


	local EfectiveKinViscosity = nil
	if self.turbViscMethod=="no" then
		EfectiveKinViscosity = KinMixViscosity
	else
		EfectiveKinViscosity = ScaleAddLinkerNumber()
		EfectiveKinViscosity:add(1.0,KinMixViscosity)
		EfectiveKinViscosity:add(Scale_RHO,KinTurbulentViscosity)

	end




	local MixViscosity = ScaleLinker();
	MixViscosity:set_import_1(Density)
	MixViscosity:set_import_2(EfectiveKinViscosity)



	local DiffusionViscosity = nil
	if self.turbViscMethod=="no" then
		DiffusionViscosity = self.rho_a * self.nu_a
	else
		DiffusionViscosity = self.rho_a * self.nu_a
	end


	local TurbulentViscosity = ScaleLinker();
	TurbulentViscosity:set_import_1(self.rho_a)
	TurbulentViscosity:set_import_2(KinTurbulentViscosity)

	
	self.TurbulentViscosity = TurbulentViscosity
	self.KinTurbulentViscosity = KinTurbulentViscosity
	self.MixViscosity = MixViscosity
	self.KinMixViscosity = KinMixViscosity
	self.EfectiveKinViscosity = EfectiveKinViscosity

	---------------------------------------------------------------------- Sediment Velocity

	local iter = 0
	local Vs = InterfaceValues:RelVel_ext(0.0,self.rho_a,self.dp,self.rho_s,self.gravity)
	local re = InterfaceValues:RE(self.nu_a*self.rho_a,self.rho_a,self.dp,Vs)
	local Cd = InterfaceValues:CD(re,self.drag_mod)


	local RelVel = RelativeVelocity(approxSpace,u)
	RelVel:set_phase_parameters(InterfaceValues)
	
	self.RelVel = RelVel

	---------------------------------------------------------------------- Diffusion (Erotion)

	local Diff_beta = 0.217
	local d1= 250e-06
	local A1= 1.673
	local k0 = 1.0 + A1 * (1-self.dp/d1)
	local Diff_factor = Diff_beta*k0*k0*Vs/(2*self.rho_a*(math.abs(self.gravity)))
	local Diffusion = GranularDiffusionLinker();
	Diffusion:set_mix_viscosity(DiffusionViscosity)
	Diffusion:set_diff_factor(Diff_factor)
	Diffusion:set_phase_parameters(InterfaceValues)

	self.Diffusion = Diffusion
	
	---------------------------------------------------------------------- Normal
	
	local ss_value = math.atan(self.FricMu_2)*180/3.1415926
	
	local Normal = DuneNormal(approxSpace,u)
	Normal:set_theta(ss_value)
	Normal:set_gradient_limit(self.grad_limit)
	Normal:set_phase_parameters(InterfaceValues)
	
	self.Normal = Normal
	
	---------------------------------------------------------------------- Vel-Diffusion (Avalanching)
	
	

	local SlipDiff = nil
	local SlipVel = nil
	if self.boolSlipDiff then
		SlipDiff = SlipDiffusion(approxSpace,u)
		SlipDiff:set_theta(ss_value)
		SlipDiff:set_diff(self.SlipVelValue)
		SlipDiff:set_slope_limit(self.slope_limit)
		SlipDiff:set_normal(Normal)
		SlipDiff:set_phase_parameters(InterfaceValues)
		
		if self.boolAveDiff then
			SlipDiff:set_diffusion(Diffusion)
		end
	else if self.boolSlipVel then
			SlipVel = SlipVelocity(approxSpace,u)
			SlipVel:set_theta(ss_value)
			SlipVel:set_vel(self.SlipVelValue)
			SlipVel:set_slope_limit(self.slope_limit)
			SlipVel:set_phase_parameters(InterfaceValues)
			SlipVel:set_normal(Normal)
		end
	end
	
	self.SlipDiff = SlipDiff
	self.SlipVel = SlipVel
	
	
	print("Clousures Discretization DONE")


end

--------------------------------------------------------------------------------
-- Discretization
--------------------------------------------------------------------------------
myProblem.Discretization = function (self,Inner_total)
	
	local NavierStokesDisc = NavierStokesFV1M (self.fct_cmp_tbl, Inner_total)
	NavierStokesDisc:set_exact_jacobian (self.bExactJac)
	NavierStokesDisc:set_stokes (self.bStokes)
	NavierStokesDisc:set_laplace (self.bNoLaplace)
	NavierStokesDisc:set_upwind (self.upwind_m)
	NavierStokesDisc:set_upwind_vol(self.upwind_t)
	NavierStokesDisc:set_peclet_blend (self.bPecletBlend)
	NavierStokesDisc:set_stabilization (self.stab, self.diffLength)
	NavierStokesDisc:set_div_correction (self.div_correction)
	NavierStokesDisc:set_transport_ip_velocity(self.boolIPVelocity)
	NavierStokesDisc:set_transport_jac(self.boolTransportJac)
	NavierStokesDisc:set_mass_term(self.boolMassTerm)
	NavierStokesDisc:set_mass_mean(self.boolDensityMean)


	NavierStokesDisc:set_density(self.Density,true)
	if self.timeMethod == "limex" and self.boolMassTerm and not(self.bStokes) then
		NavierStokesDisc:set_limex_correction(true)
	end

		
	if self.boolRelativeVel then
		NavierStokesDisc:set_relative_velocity(self.RelVel,self.riemman)
		NavierStokesDisc:set_upwind_rel(self.upwind_r)
	end
	if self.boolSlipDiff then
		NavierStokesDisc:set_diffusion(self.SlipDiff)
	else
		if self.boolAveDiff then
			NavierStokesDisc:set_diffusion(self.Diffusion)
		end
		if self.boolSlipVel then
			NavierStokesDisc:set_slip_velocity(self.SlipVel)
		end
		
	end
	if(self.boolpress_jump) then
		NavierStokesDisc:set_pressure_jump ( self.diffLength)
	end

	NavierStokesDisc:set_kinematic_viscosity (self.EfectiveKinViscosity)

	NavierStokesDisc:set_average_gamma(self.gamma)
	NavierStokesDisc:set_phase_parameters(self.InterfaceValues)
	
	print("Space Discretization DONE")
	
	self.NavierStokesDisc = NavierStokesDisc
	return NavierStokesDisc
end

--------------------------------------------------------------------------------
-- Time Discretization
--------------------------------------------------------------------------------
myProblem.TimeDiscretization = function (self,domainDisc)
	
	print("Time Discretization")
	-- create the assembled operator for the solver
	local timeDisc = nil
	-- create time discretization
	if self.timeMethod=="cn" then
		timeDisc = ThetaTimeStep(domainDisc)
		timeDisc:set_theta(0.5) -- Crank-Nicolson method
	end
	if self.timeMethod=="euler" then
		timeDisc = ThetaTimeStep(domainDisc)
		timeDisc:set_theta(1) -- implicit Euler
	end
	if self.timeMethod=="fracstep" then
		timeDisc = ThetaTimeStep(domainDisc,"FracStep")
	end
	if self.timeMethod=="alex" then
		timeDisc = ThetaTimeStep(domainDisc, "Alexander")
	end
	
	print("Time Discretization DONE")
	
	self.timeDisc = timeDisc
	return timeDisc
end

--------------------------------------------------------------------------------
-- SOLVER
--------------------------------------------------------------------------------
myProblem.CreateSolver = function (self, domainDisc, approxSpace)
	
	--------------
	-- LineSearch
	--------------
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
	----------------------
	-- NoLinear COnvCheck
	----------------------
	local NewtonSteadyConvCheck=ConvCheck(self.max_newton_steps_steady_state, self.SteadyAbsDefect, self.SteadyRedDefect, true)
	local NewtonConvCheck=ConvCheck(self.max_newton_steps_transient, self.AbsDefect, self.RedDefect, true)
	local LinearConvCheckImp=ConvCheck(self.max_linear_steps_Imp, self.LinAbsDefectImp, self.LinRedDefectImp, true)
	local LinearConvCheckLim=ConvCheck(self.max_linear_steps_Lim, self.LinAbsDefectLim, self.LinRedDefectLim, true)
	local LimexConvCheck=ConvCheck(1, 1e-12, 1e-12, true)
	      LimexConvCheck:set_supress_unsuccessful(true)
	
	--------------
	-- Smoothers
	--------------
	-- base solver
	baseSolver = LU()
	baseSolver = AgglomeratingSolver(SuperLU());
	
	
	ilu = ILU()
	ilu:set_beta(self.value_beta)
	ilu:set_damp(self.damping_mg)
	ilu:set_ordering_algorithm(TopologicalOrdering())
	ilu:set_sort(true)
	--ilu:set_sort_eps(1.e-50)
	ilu:set_inversion_eps(1.e-8)
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
	

	------------------
	-- preconditioners
	-------------------

	gmg = GeometricMultiGrid(approxSpace)
	gmg:set_discretization(domainDisc)
	gmg:set_base_level(self.numPreRefs)
	gmg:set_base_solver(baseSolver)
	gmg:set_smoother(ilu)
	gmg:set_cycle_type(1)
	gmg:set_num_presmooth(3)
	gmg:set_num_postsmooth(3)
	gmg:set_rap( true)
	gmg:set_smooth_on_surface_rim(false)

	-- gmg:set_damp(MinimalResiduumDamping())
	-- gmg:set_damp(0.8)
	-- gmg:set_damp(MinimalEnergyDamping())
	
	
	-----------------
	-- Linear Solver
	-----------------
	
		-- create Linear Solver
	--GMresSolver = GMRES(20)
	--GMresSolver:set_preconditioner(gmg)
	--GMresSolver:set_convergence_check(LinearConvCheck)
	
	-- create Linear Solver
	BiCGStabSolverImp = BiCGStab()
	BiCGStabSolverImp:set_preconditioner(gmg)
	BiCGStabSolverImp:set_convergence_check(LinearConvCheckImp)
	
	BiCGStabSolverLim = BiCGStab()
	BiCGStabSolverLim:set_preconditioner(gmg)
	BiCGStabSolverLim:set_convergence_check(LinearConvCheckLim)

	--gmgSolver = LinearSolver()
	--gmgSolver:set_preconditioner(gmg)
	--gmgSolver:set_convergence_check(LinearConvCheck)
	
	--ilutSolver = LinearSolver()
	--ilutSolver:set_preconditioner(ilu)
	--ilutSolver:set_convergence_check(LinearConvCheck)
	
	
	-- choose a solver
	LinearSolverLim = BiCGStabSolverLim
	LinearSolverImp = BiCGStabSolverImp
	--LinearSolver = GMresSolver
	--LinearSolver = gmgSolver
	--LinearSolver = ilutSolver
	
	
	
	

	local NewtonSolverSteady = nil
	if self.doSteadyState then
		NewtonSolverSteady = NewtonSolver()
		NewtonSolverSteady:set_linear_solver(LinearSolverImp)
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
		NLSolver:set_linear_solver(LinearSolverLim)
		NLSolver:set_convergence_check(LimexConvCheck)
		NLSolver:auto_update(false)
		limex = myProblem:LimexObject( domainDisc, NLSolver)
	else
		NLSolver:set_linear_solver(LinearSolverImp)
		NLSolver:set_convergence_check(NewtonConvCheck)
		NLSolver:set_line_search(NewtonLineSearch)
		NLSolver:set_reassemble_J_freq(0)
		NLSolver:auto_update(true)
		op = AssembledOperator(self.timeDisc)
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
		TransientNewtonUpdater:set_min_threshold(-1.1)
		NLSolver:setNewtonUpdater(TransientNewtonUpdater)
	end
	
	
	self.NewtonSolverDescSteady = NewtonSolverDescSteady
	self.NewtonSolverDesc = NewtonSolverDesc
	self.boolSolverDesc = true
	
	print("Solver Setting DONE")
	return op, NLSolver, NewtonSolverSteady, limex
end


--------------------------------------------------------------------------------
-- OutputParameters
--------------------------------------------------------------------------------
myProblem.OutputParameters = function (self)

	local out = VTKOutput()
	out:clear_selection()
	out:select_all(false)
	if self.dim == 2 then
		out:select_nodal ("u,v", "velocity")
	else
		out:select_nodal ("u,v,w", "velocity")
	end
	out:select_nodal ("u", "u")
	out:select_nodal ("v", "v")
	if self.dim == 3 then
		out:select_nodal ("w", "w")
	end
	out:select_nodal ("p", "p")
	out:select_nodal ("c", "c")
	out:select(self.Density, "Rho")
	out:select(self.NavierStokesDisc:einstein_viscosity(), "Mu_eins")
	out:select(self.NavierStokesDisc:mix_viscosity(), "Mu_I")
	out:select(self.MixViscosity, "MixViscosity")
	out:select(self.TurbulentViscosity, "Mu_turb")
	
	out:select(self.RelVel, "RelVel")
	out:select(self.NavierStokesDisc:particle_pressure(), "Ps")
	out:select(self.NavierStokesDisc:particle_pressure_grad(), "DPs")
	out:select(self.NavierStokesDisc:velocity_grad(), "Gamma")
	out:select(self.gamma, "MeanGamma")
	if (self.boolSlipDiff) then
		out:select_element(self.SlipDiff, "SDiff")
	else
		if self.boolSlipVel then
			out:select_element(self.SlipVel, "SVel")
		end
	end
	out:select_element(self.Diffusion, "D")
	out:select_element(self.Normal, "n")

	print("Output file setting DONE")
	return out
end
--------------------------------------------------------------------------------
-- CheckPoint
--------------------------------------------------------------------------------
myProblem.SaveCheckPoint = function (self,u,folder)

	SaveToFile(u, folder .. "/CheckPoint" .. ".vec")
	print("Saving CheckPoint: DONE")

end

myProblem.LoadCheckPoint = function (self,u,folder)

	
	local boolInterpolate = false
	local filename = folder .. "/Integral.txt"

	local step = nil
	local file = io.open(filename, "r")
	local time_work_total = 0.0

	if file then
		print("Loading CheckPoint")
		for line in file:lines() do

			-- Skip the header
			if not line:match("^Step") and not line:match("^%-") then

				local columns = {}

				for value in line:gmatch("%S+") do
					table.insert(columns, value)
				end

				-- First column = Step
				step = tonumber(columns[1])

				-- Sixth column = TNSteps
				time_work_total = tonumber(columns[6])
			end
		end

		file:close()
		
		ReadFromFile(u, folder .. "/CheckPoint" .. ".vec")
	else
		boolInterpolate = true
		step = 0
	end

	
	time = self.DT * step
	
	
 return time, step, time_work_total, boolInterpolate
end
--------------------------------------------------------------------------------
-- Writing Output parameters
--------------------------------------------------------------------------------

myProblem.WriteValues = function (self, folder, step, time, Value_inner1, Value_inner2, WorkTime, TotalWorkTime, Newton_Steps, Newton_Steps_fail,linsolver_calls,linsolver_steps,boolTotal)

	if(boolTotal) then
		file = io.open(folder .. "/Integral.txt", "a")
		file:write(string.format("-----------------------------------------------------------------------------------------------------------\n"))
		file:close()
	end
	if(step == 0) then
		file = io.open(folder .. "/Integral.txt", "w+")
		file:write(string.format("Step\tTime\t\tVol-Dom_1\tVol-Dom_2\tWork time\tTotal work time\tTNSteps\tSNSteps\tFNSteps\tLinCalls LinSteps\n"))
		file:write(string.format("%d\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%d\t%d\t%d\t%d\t%d\n", step, time, Value_inner1, Value_inner2, WorkTime, TotalWorkTime, Newton_Steps, Newton_Steps-Newton_Steps_fail, Newton_Steps_fail,linsolver_calls,linsolver_steps))
		file:close()
	else
		file = io.open(folder .. "/Integral.txt", "a")
		file:write(string.format("%d\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%d\t%d\t%d\t%d\t%d\n", step, time, Value_inner1, Value_inner2, WorkTime, TotalWorkTime, Newton_Steps, Newton_Steps-Newton_Steps_fail, Newton_Steps_fail,linsolver_calls,linsolver_steps))
		file:close()
	end
	if(boolTotal) then
		file = io.open(folder .. "/Integral.txt", "a")
		file:write(string.format("-----------------------------------------------------------------------------------------------------------\n"))
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
	local tBefore_s= os.clock()
	if not solver:apply(u) then
		print("===> THE PREPARATION PHASE FAILED! <===")
	end
	local tAfter_s = os.clock()
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
myProblem.SolveNonlinearProblem = function (self, u, solver, op, solTimeSeries, dt, step, StartTime, EndTime)

	                    
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
		self.timeDisc:prepare_step(solTimeSeries, do_dt)
	
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
					self.timeDisc:prepare_step(solTimeSeries, do_dt)
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
            
                time= self.timeDisc:future_time()                                        -- update new time
                                                                                
                oldestSol = solTimeSeries:oldest()                                  -- get oldest solution
                                                                                    
				VecScaleAssign(oldestSol, 1.0, u)                                             -- copy values into oldest solution (we reuse the memory here)
                                                                                    
                solTimeSeries:push_discard_oldest(oldestSol, time)                  -- push oldest solutions with new values to front, oldest sol pointer is poped from end
                
                CompletedStep = true
                
                
            else
                
                time2 = self.timeDisc:future_time()                                      -- update new time
                                                                                    
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
	local dtlimex = self.DTLimex

	--  Euclidean (algebraic) norm
	--local estimator = Norm2Estimator()



	--print (estimator)
	local limexEstimator = CompositeGridFunctionEstimator()
	
	if self.VelErrorNorm == "L2" then
		limexEstimator:add(L2ComponentSpace("u", 2, self.alphaVel))
		limexEstimator:add(L2ComponentSpace("v", 2, self.alphaVel))
	elseif self.VelErrorNorm == "H1" then
		limexEstimator:add(H1SemiComponentSpace("u", 2 , ConstUserMatrix(self.alphaVel)))
		limexEstimator:add(H1SemiComponentSpace("v", 2 , ConstUserMatrix(self.alphaVel)))
	else
		print ("LimexErrorEstimator for velocity not defined"); exit();
	end
	
	if self.PressErrorNorm == "L2" then
		limexEstimator:add(L2ComponentSpace("p", 2, self.alphaPress))
	elseif self.PressErrorNorm == "H1" then
		limexEstimator:add(H1SemiComponentSpace("p", 2, ConstUserMatrix(self.alphaPress) ))
	else
		print ("LimexErrorEstimator for Pressure not defined"); exit();
	end
	
	if self.VolErrorNorm == "L2" then
		limexEstimator:add(L2ComponentSpace("c", 2, self.alphaVol))
	elseif self.VolErrorNorm == "H1" then
		limexEstimator:add(H1SemiComponentSpace("c", 2, self.alphaVol))
	else
		print ("LimexErrorEstimator for VolumeFraction not defined"); exit();
	end
	

	
	limexEstimator:use_strict_relative_norms(1)
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
	  dtred = self.red_factor_fail,
	  dtIncr = self.incr_factor,
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
myProblem.SolveNonlinearProblemLimex = function (self, u, limex, NLSolver, time_step, StartTime, EndTime, NewtonLimexSteps)
	
	if(time_step == 1) then
		limex:set_dt_min(1e-012)
	end
	limex:set_start_step(1)
	limex:apply(u, EndTime, u, StartTime)
	n_step = limex:get_step()-1
	
    local Newton_Steps = NLSolver:total_linsolver_calls()/(self.nstages+1)
	local Newton_Steps_fail = 0/(self.nstages+1)
	linsolver_calls_step = NLSolver:total_linsolver_calls()
	linsolver_steps_step =  NLSolver:total_linsolver_steps()
	limex:set_time_step(self.DTmax/(math.max(n_step-1,1)))
	
	if(time_step == 1) then
		limex:set_dt_min(self.DTmin)
	end
	
	NLSolver:clear_average_convergence();

  return Newton_Steps, Newton_Steps_fail, linsolver_calls_step, linsolver_steps_step, 1
end


--------------------------------------------------------------------------------
-- ParaView Contour Extraction
--------------------------------------------------------------------------------
myProblem.RunParaViewContour = function(self, rank, folder_vtk)

    print(
        "RunParaViewContour ENTER: rank=" .. tostring(rank) ..
        ", boolData=" .. tostring(self.boolData) ..
        ", folder=" .. tostring(folder_vtk)
    )

    if not self.boolData then
        print("RunParaViewContour SKIPPED: boolData=false")
        return
    end

    if rank ~= 0 then
        return
    end

    print("============================================")
    print("Starting ParaView contour extraction")
    print("============================================")

    local command =
        "pvpython FunctionTools.py contour " ..
        "\"" .. folder_vtk .. "\" " ..
        tostring(self.dim) .. " " ..
        "\"" .. self.data_name .. "\""

    print("Input folder : " .. folder_vtk)
    print("Dimension    : " .. tostring(self.dim))
    print("Data name    : " .. self.data_name)
    print("Command      : " .. command)

    local result = os.execute(command)

    print("ParaView return code:", result)

    if result ~= 0 then
        print("WARNING: ParaView contour extraction failed.")
    else
        print("ParaView contour extraction completed successfully.")
    end

    print("============================================")
    print("ParaView contour extraction finished")
    print("============================================")
end

return myProblem


