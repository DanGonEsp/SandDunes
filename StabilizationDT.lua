------------------------------------------------------------------------------------------
-- Navier-Stokes equation, 3d
-- Discretization: Vertex-centered, stabilized
------------------------------------------------------------------------------------------

-- Load utility scripts (e.g. from from ugcore/scripts)
ug_load_script ("ug_util.lua")
ug_load_script ("util/load_balancing_util.lua")

ug_load_script("util/domain_disc_util.lua")
ug_load_script("navier_stokes_util.lua")
ug_load_script("util/conv_rates_kinetic.lua")


------------------------------------------------------------------------------------------
-- Get command line parameters
------------------------------------------------------------------------------------------

-- Physical parameters
dim 		= util.GetParamNumber("-dim", 2, "dimensionality of the problem")
geometry	= util.GetParam ("-geom", "Dune2D_tri_5")
nu_a 	= util.GetParamNumber("-visc_a", 1.48e-05, "kinematic viscosity")
nu_s 	= util.GetParamNumber("-visc_s", 1.48e-03, "kinematic viscosity")
rho_a 	= util.GetParamNumber("-rho_a", 1.2, "Air Density")
rho_s 	= util.GetParamNumber("-rho_s", 1000, "Sand Density")
dp 	= util.GetParamNumber("-diameter", 0.5e-03, "Particle Diameter")
inflow		= util.GetParamNumber("-inflow", 1.0, "max. inflow velocity")
c_init		= util.GetParamNumber("-initial concentration", 1.0, "max volume fraction")


alpha_max		= util.GetParamNumber("-max_concentration", 0.6, "max volume fraction")
packing_factor		= util.GetParamNumber("-packing_facotr", 0.5999, "max volume fraction")
alpha_min		= util.GetParamNumber("-min concentration", 0.5, "max volume fraction")
viscosity_model		= util.GetParamNumber("-granular_model", 2, "Options:  0 Constant, 1 Proportional, 2 Einstein model, 3 Rheology(I), 4 Rheology(I) + Einstein model")
density_model  = util.GetParam("-density_model", "linear", "constant, linear")
interface_value  = util.GetParamNumber("-interface_value", 0.5, "interface value")


-- Numerical parameters
harmonic=false
interface_harmonic=false

jumpPressure = false



dt = util.GetParamNumber("-dt",0.05)
numTimeSteps =  util.GetParamNumber("-numTimeSteps", 20	)
EndTime = util.GetParamNumber("-EndTime", 6.2, "EndTime")
boolEndTime 	= util.GetParamNumber("-boolEndTime", true)
outputFactor     = util.GetParam("-output", 1, "output every ... steps")
max_newton_steps_1=util.GetParamNumber("-numNewtonSteps", 100)
max_newton_steps_2=util.GetParamNumber("-numNewtonSteps", 100)
max_linear_steps=util.GetParamNumber("-numLinearIter", 200)
StatBool = util.GetParamNumber("-StatBool", false, "Stationary state")


numPreRefs = util.GetParamNumber("-numPreRefs", 2)
numRefs = util.GetParamNumber("-numRefs",2)
timeMethod = util.GetParam("-timeMethod","euler")

CFL_max= util.GetParamNumber("-cfl", 100, "max  CFL number")
DT_max= util.GetParamNumber("-DT_max", 2, "max  DT")
DT_min= util.GetParamNumber("-DT_min", 0.001, "min  DT")
UpdateDt 	= util.GetParamNumber("-UpdateDT", 5)
modifyDT 	= util.GetParamNumber("-modifyDT", false)
incr_factor 	= util.GetParamNumber("-CFL_factor", 1.1)
red_factor 	= util.GetParamNumber("-CFL_factor", 0.9)



turbViscMethod = util.GetParam("-turbulenceModel","no")
modellconstant = util.GetParamNumber("-c",0.1)

bStokes 	= util.GetParam("-Stokes", false ,"If defined, only Stokes Eq. computed")
bNoLaplace 	= util.GetParam("-noLaplace", false,"If defined, only laplace term used")
bExactJac 	= util.GetParam("-exactJac", false,"If defined, exact jacobian used")
bPecletBlend= util.GetParam("-PecletBlend", false,"If defined, Peclet Blend used")
upwind      = util.GetParam("-upwind", "full", "Upwind type full or lps")
diffLength  = util.GetParam("-difflength", "cor", "fivepoint, raw, corDiffusion length type")
bPac        = util.GetParam("-pac", false,"If defined, pac upwind used")
if (jumpPressure) then
	file_name ="Pressure2"
	stab        = util.GetParam("-stab", "fields_2", "Stabilization type (fields or flow viscosity or karimian)")
else
	file_name ="NoPressure"
	stab        = util.GetParam("-stab", "fields_2", "Stabilization type (fields or flow viscosity or karimian)")
end



vtk_file_name = file_name .. "-lev" .. numRefs
if bStokes then
	vtk_file_name = vtk_file_name .. "-Stokes"
end


-- Grid file name
gridName = geometry .. ".ugx"

-- Subsets used in the problem
allSubsets = "Inner,Left, Right,Top, Bottom,"

print (" Geometry: " .. geometry .. " (file " .. gridName .. "), dim = " .. dim)
print (" Physical parameter:")
print ("	inflow		= " .. inflow)
print ("	Stokes		= " .. tostring (bStokes))
print (" Numerical parameter:")
print ("	numRefs		= " .. numRefs)
print ("	numPreRefs	= " .. numPreRefs)
print ("	noLaplace		= " .. tostring (bNoLaplace))
print ("	exactJac	= " .. tostring (bExactJac))
print ("	PecletBlend 	= " .. tostring (bPecletBlend))
print ("	upwind	= " .. upwind)
print ("	pac			= " .. tostring (bPac))
print ("	stab		= " .. stab)
print ("	difflength	= " .. diffLength)

------------------------------------------------------------------------------------------
-- Initialize UG4, load, refine and distribute the grid
------------------------------------------------------------------------------------------

InitUG (dim, AlgebraType("CPU", dim +2))

if dim == 3 then
	fct_cmp_tbl = {"u", "v", "w", "p"}
	vel_cmp_tbl = {"u", "v", "w"}
else
	fct_cmp_tbl = {"u", "v", "p"}
	vel_cmp_tbl = {"u", "v"}
end	

-- Create the domain, load the grid and refine it
dom = util.CreateDomain (gridName, numPreRefs)
balancer.RefineAndRebalanceDomain (dom, numRefs - numPreRefs)

print ("Domain-info:")
print (dom:domain_info():to_string())

-- Create the vertex-centered approximation space
approxSpace = ApproximationSpace (dom)

approxSpace:add_fct("u", "Lagrange",1,allSubsets)
approxSpace:add_fct("v", "Lagrange",1,allSubsets)
--approxSpace:add_fct("w", "Lagrange",1,allSubsets)
approxSpace:add_fct("p", "Lagrange",1,allSubsets)
approxSpace:add_fct("c", "Lagrange",1,allSubsets)


approxSpace:init_levels()
approxSpace:init_top_surface()
approxSpace:print_statistic()

util.solver.defaults.approxSpace = approxSpace


-- Order the DoFs:
OrderLex (approxSpace,  "x")
--OrderCuthillMcKee(approxSpace,true)
-- grid function for the solution
u = GridFunction(approxSpace)
u:set(0)

--------------------------------
--------------------------------
-- Lua Functions
--------------------------------
--------------------------------
H_domain=5
H_0=2
h_0=0.3
sigma=3--0.4
mu_c=15
L2=5
L=6--2.5
--[[
H_0=1.5
h_0=0.1
sigma=5
mu_c=14
L2=0*4	
L=10--2.5]]


mu_c1=mu_c-L2/2
mu_c2=mu_c+L2/2
sigma1=sigma
sigma2=sigma
d=1.2*math.pow(1/2,numRefs-1)
ss=1.0
k1=0.05 
k2=0.02 
---------------------------------------------------------------------- Initial DuneShape
function Dune(x) 
	if (x< mu_c1) then
		return H_0*   math.exp(     -0.5*math.pow((x-mu_c1)/sigma1,2) ) -h_0+(ss-math.exp(-k1*(x-mu_c1)))
	else
		if (x< mu_c2) then
			return H_0-h_0
		else
			return H_0*   math.exp(     -0.5*math.pow((x-mu_c2)/sigma2,2) ) -h_0  -(ss-math.exp(-k2*(x-mu_c2)))
			
		end
	end
end
function Dune1(x) 
	if (x< mu_c1) then
		return -(Dune(x)+h_0)*(x-mu_c1)/math.pow(sigma1,2)+(k1*math.exp(-k1*(x-mu_c1)))
	else
		if (x< mu_c2) then
			return 0
		else
			return -(Dune(x)+h_0)*(x-mu_c2)/math.pow(sigma2,2)+(-k2*math.exp(-k2*(x-mu_c2)))			
		end
	end
end
function Dune2(x) 
	if (x< mu_c1) then
		return (Dune(x)+h_0)*(math.pow((x-mu_c1)/math.pow(sigma1,2),2)-1/math.pow(sigma1,2))-(k1*k1*math.exp(-k1*(x-mu_c1)))
	else
		if (x< mu_c2) then
			return 0
		else
			return  (Dune(x)+h_0)*(math.pow((x-mu_c2)/math.pow(sigma2,2),2)-1/math.pow(sigma2,2))+(k2*k2*math.exp(-k2*(x-mu_c2)))			
		end
	end
end

function dist(x0,y0)
	e=1
	x1=20
	x2=20
	y1=Dune(x0)
	n=1
	while (e>1e-04 and n<100) do
		x1=x2
		y1=Dune(x1)
		x2=x1-(x0-x1  + (y0-y1)*Dune1(x1))/(Dune2(x1)*(y0-Dune(x1))-(math.pow(Dune1(x1),2)+1) )

		e=math.abs(x2-x1)
		n=n+1
	end
	return math.sqrt(math.pow(x2-x0,2)+math.pow(Dune(x2)-y0,2))
end

dx=math.pow(1/2,(numRefs+2))

	
---------------------------------------------------------------------- Initial Velocity

--[[function StartValueX(x,y) 
	interface=0
	dd=dist(x,y)
	ds=dx
	--ds=0.1*dx
	if dd<ds then
		interface=1
	end
	if (y>Dune(x,y) and y>0.5*dx) then
		if interface==0 then
			hh=14.1856
			nn=2.5
			cc=1.0--math.pow(y/hh,5)
			return inflow*(math.min(1.0, math.pow(y/hh,1/nn))*(1-cc) +(cc)* (hh - y) * (4*y ) / (hh * hh))
		else
			return 0.0
		end
	else
		return 0.0
	end
end]]
function StartValueX(x,y) 
	hh=14.1856
	nn=2.5
	cc=1.0--math.pow(y/hh,5)
	return inflow*(math.min(1.0, math.pow(y/hh,1/nn))*(1-cc) +(cc)* (2*hh - y) * (y ) / (hh * hh))
end
function StartValueY(x,y) 
	return 0
end
function StartValueZ(x,y) 
	return 0
end
---------------------------------------------------------------------- Initial Pressure
Pstd=0.0	
function StartValueP(x,y) 
	if y>Dune(x) then
		return Pstd
	else 
	return  0*9.81*rho_s*c_init*(Dune(x)-y)+Pstd
	end
end

---------------------------------------------------------------------- Initial VolumeFraction 



function VolumeFraction(x,y)
	interface_o=0
	interface_i=0
	dd=dist(x,y)
	ds_o=4*dx
	ds_i=4*dx
	if (dd<ds_o) then
		interface_o=1
	end
	if dd<ds_i then
		interface_i=1
	end
	if y>Dune(x,y) then
		if interface_o==0 then
			return 0.00
		else
			return interface_value*(1.0 - dd/ds_o)
			--return 0.0
		end
	else
		if interface_i==0 then
			return  c_init
		else
			return interface_value*(1.0 - dd/ds_o)+ c_init*dd/ds_i
			--return math.min(1.0, interface*(dd/(ds-dd)+1.0))
			--return c_init 
		end
	end
	
	
end
function VolumeFraction2(x,y)
	dd=dist(x,y)
	ds=20*dx
	k=1.0/dx
	
	if y>Dune(x,y) then
		if dd<ds then
			dd=-dd
			k=k
			return c_init / (1.0 + math.exp(-k*dd))
		else
			return 0.0
		end
	else
		if dd<ds then
			k=k
			return c_init / (1.0 + math.exp(-k*dd))
		else
			return c_init 
		end
	end
	
	
	
end
function InitialValue_FractionVolume(x,y)
	--if y>0 then 
		value= VolumeFraction(x,y) 
	--else
		--value = BoundaryVolumeFraction(x,y)
	--end
	return value
end
---------------------------------------------------------------------- Boundary Condition  

----------------------------------------------------------- Bottom

function BoundaryVolumeFraction(x,y)
	if x>2 and x<28 then--y>Dune(x,y) then
		return c_init
	else
		return 0.0
	end
end

function BottomFlux(x,y) return 0 end

----------------------------------------------------------- Outlet

function PressureOutlet(x,y,t)
	return 0*rho_a*11.83*(H_domain-y)
end

function NeuBND(c) return inflow*c end
function DNeuBND(c) return inflow end
Flux = LuaUserFunctionNumber("NeuBND", 1);
Flux:set_deriv(0, "DNeuBND")


---------------------------------------------------------------------- Inflow
function inflowVel3d_2(x, y, t)
	hh=0.1
	return inflow*math.min(1.0, math.pow(y/hh,4)), 0.0

end
---------------------------------------------------------------------- Inflow
function inflowVel3d(x, y, t)
	return StartValueX(x,y), StartValueY(x,y)
end
---------------------------------------------------------------------- Density


Density = GranularDensityLinker(); 
Density:set_fluid_density(rho_a)
Density:set_particle_density(rho_s*packing_factor)
Density:set_harmonic_value(harmonic)
Density:set_surface_value(interface_harmonic)
Density:set_model(density_model)
Density:set_interface_volume_fraction(interface_value)
---------------------------------------------------------------------- Viscosity

Ps = ParticlePressureLinker(); 
Ps:set_particle_diameter(dp)
Ps:set_particle_density(rho_s)
Ps:set_fluid_Visc(nu_a*rho_a)
Ps:set_alpha_max(alpha_max)
Ps:set_alpha_min(alpha_min)
Ps:set_packing_factor(packing_factor)



Visc = GranularViscosityLinker(); 
Visc:set_granular_model(viscosity_model)
Visc:set_particle_diameter(dp)
Visc:set_particle_density(rho_s)
Visc:set_particle_kinematicVisc(nu_s)
Visc:set_fluid_density(rho_a)
Visc:set_fluid_Visc(nu_a*rho_a)
Visc:set_alpha_max(alpha_max)
Visc:set_alpha_min(alpha_min)
Visc:set_packing_factor(packing_factor)
Visc:set_mix_density(Density)
Visc:set_harmonic_value(harmonic)
Visc:set_surface_value(interface_harmonic)
Visc:set_interface_volume_fraction(interface_value)

PJump=PressureJumpLinker()
PJump:set_mix_density(Density)
PJump:set_kinematic_viscosity(Visc)

PjumpShape= JumpShapeLinker()
PjumpShape:set_interface_volume_fraction(interface_value)


normal = InterfaceNormalLinker()
normal:set_interface_volume_fraction(interface_value)


---------------------------------------------------------------------- Gravity

Gravity = ConstUserVector(0.0)
Gravity:set_entry(1, -9.81)


------------------------------------------------------------------------------------------
-- Compose the discretization
------------------------------------------------------------------------------------------

-- inner space
NavierStokesDisc = NavierStokesFV1(fct_cmp_tbl, {"Inner"})

NavierStokesDisc:set_peclet_blend(bPecletBlend) 
NavierStokesDisc:set_exact_jacobian(bExactJac)
NavierStokesDisc:set_laplace(bNoLaplace)
NavierStokesDisc:set_stokes(bStokes)
NavierStokesDisc:set_upwind (upwind)
NavierStokesDisc:set_stabilization (stab, diffLength)
NavierStokesDisc:set_pac_upwind (bPac)



---------------------------------------------------------------------------------------
-- Set properties input
---------------------------------------------------------------------------------------
NavierStokesDisc:set_density(Density)
NavierStokesDisc:set_kinematic_viscosity(Visc)
NavierStokesDisc:set_density_ref(rho_a)

if (jumpPressure) then
	VolFraction = GridFunctionNumberData(u, "c");
	NavierStokesDisc:set_jump_shape(PjumpShape,diffLength)
	NavierStokesDisc:set_interface_normal(normal)
	NavierStokesDisc:set_vol_fraction(VolFraction)
	NavierStokesDisc:set_interface_value(interface_value)
end

--NavierStokesDisc:set_source(Gravity)
--NavierStokesDisc:set_source(RhoGrad)
--NavierStokesDisc:set_source_surface(Fs)
--NavierStokesDisc:set_relative_velocity(W)





InletDisc = NavierStokesInflow (NavierStokesDisc)
InletDisc:add ("inflowVel3d", "Left,Top")

-- boundary condition at the outlet
OutletDisc = NavierStokesNoNormalStressOutflow (NavierStokesDisc)
OutletDisc:add ("Right")

-- boundary condition at the impermeable walls
WallDisc = NavierStokesWall (NavierStokesDisc)
WallDisc:add (" Bottom")
print("Navier Stokes Equation created.")

velocity = GridFunctionVectorData(u, "u,v");


TransportEq = ConvectionDiffusion("c", "Inner", "fv1")
TransportEq:set_velocity(NavierStokesDisc:velocity_ip())
--TransportEq:set_velocity(NavierStokesDisc:velocity())
TransportEq:set_diffusion(0)
TransportEq:set_upwind(UpwindFV1(upwind)) 


-- create dirichlet boundary for concentration
dirichletBND = DirichletBoundary()
dirichletBND:add("BoundaryVolumeFraction", "c", "Bottom")
--dirichletBND:add(0.001, "c", "Bottom")

OutflowBND = ConvectionDiffusionOutflowFV1(TransportEq)
OutflowBND:add( "Right, Top, Left, Bottom")	

NeumannBND = NeumannBoundaryFV1("c")
--NeumannBND:add( BottomFlux,"Bottom,Top,Left","Inner")	
--NeumannBND:add( Flux,"Right","Inner")


print("Transport Equation created.")


---------------------------------------------------------------------------------------
-- Viscosity, Density and Gravitation Input
---------------------------------------------------------------------------------------


--Visc:set_input(0, TransportEq:value())

Visc:set_volume_fraction(TransportEq:value())
Visc:set_velocity_gradient(NavierStokesDisc:velocity_grad())
Visc:set_particle_pressure(NavierStokesDisc:pressure())

Density:set_volume_fraction(TransportEq:value())

normal:set_volume_fraction(TransportEq:value())
normal:set_volume_grad(TransportEq:gradient())


PJump:set_volume_fraction(TransportEq:value())
PJump:set_volume_grad(TransportEq:gradient())
PJump:set_velocity_gradient(NavierStokesDisc:velocity_grad())
PJump:set_pressure_grad(NavierStokesDisc:pressure_grad())



PjumpShape:set_volume_fraction(TransportEq:value())




domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(NavierStokesDisc)
domainDisc:add(WallDisc)
domainDisc:add(InletDisc)
domainDisc:add(OutletDisc)


domainDisc:add(TransportEq)
--domainDisc:add(NeumannBND)
--domainDisc:add(dirichletBND)
--domainDisc:add(OutflowBND)

-- create operator from discretization





-- create the assembled operator for the solver

-- create time discretization
if timeMethod=="cn" then
	timeDisc = ThetaTimeStep(domainDisc)
	timeDisc:set_theta(0.5) -- Crank-Nicolson method
end
if timeMethod=="euler" then
	timeDisc = ThetaTimeStep(domainDisc)
	timeDisc:set_theta(1) -- implicit Euler
end	
if timeMethod=="fracstep" then
	timeDisc = ThetaTimeStep(domainDisc,"FracStep")
end
if timeMethod=="alex" then
	timeDisc = ThetaTimeStep(domainDisc, "Alexander")
end

op = AssembledOperator(timeDisc)


op:init()

------------------------------------------------------------------------------------------
-- Set up the solver
------------------------------------------------------------------------------------------

-- For debugging only (to write the intermediate data): --
--util.debug_dir = "FLOW_DEBUG"
--util.debug = { vtk = true, conn_viewer = true }
--util.CreateGridFuncDebugWriter (approxSpace)
-- --



-- create algebraic Preconditioner
ilusmoother = ILU()
ilusmoother:set_beta(1e-4)
--ilusmoother:set_inversion_eps(1e-30)
--ilu:set_debug(dbgWriter)

ilutsmoother = ILUT()
ilutsmoother:set_threshold(1e-4)
--ilutsmoother:set_info(true)

local smooth = util.smooth.parseParams()
smoother_local = util.smooth.create(smooth)

ilutSolver = LinearSolver()
ilutSolver:set_preconditioner(ilusmoother)
ilutSolver:set_convergence_check(ConvCheck(10000, 1e-6, 1e-1, false))

exactSolver = LU()

baseSolver=exactSolver
--baseSolver=ilutSolver


--smoother=vanka
--smoother=smoother_local
smoother=ilusmoother
smoother:set_damp(1)

numPreSmooth=1
numPostSmooth=1
gmg = GeometricMultiGrid(approxSpace)
gmg:set_discretization(domainDisc)
gmg:set_base_level(numPreRefs)
gmg:set_base_solver(baseSolver)
gmg:set_smoother(smoother)
gmg:set_cycle_type("V")
gmg:set_num_presmooth(numPreSmooth)
gmg:set_num_postsmooth(numPostSmooth)
gmg:set_damp(MinimalResiduumDamping())
--gmg:set_damp(0.1)
gmg:set_rap(true)
--gmg:set_ordering("x")
--gmg:set_damp(MinimalEnergyDamping())

--gmg:set_debug(dbgWriter)

-- create Linear Solver
convCheck = ConvCheck(max_linear_steps, 1e-15, 1e-14, true)
convCheck:set_verbose(true)


BiCGStabSolver = BiCGStab()
BiCGStabSolver:set_preconditioner(gmg)
--BiCGStabSolver:set_preconditioner(vanka)
BiCGStabSolver:set_convergence_check(convCheck)
BiCGStabSolver:set_compute_fresh_defect_when_finished(true)




gmgSolver = LinearSolver()
gmgSolver:set_preconditioner(gmg)
gmgSolver:set_convergence_check(convCheck)
-- create Exact solver


-- choose a solver
solver = exactSolver
solver = baseSolver
solver = BiCGStabSolver
--solver = vankaSolver
--solver = gmgSolver --- este era el bueno
--solver = ilutSolver
--solver = linSolver


newtonConvCheck = ConvCheck()
newtonConvCheck:set_maximum_steps(max_newton_steps_1)
newtonConvCheck:set_minimum_defect(1e-04)
newtonConvCheck:set_reduction(1e-10)
newtonConvCheck:set_verbose(true)

newtonLineSearch = StandardLineSearch()
newtonLineSearch:set_maximum_steps(5)
newtonLineSearch:set_lambda_start(1)
newtonLineSearch:set_reduce_factor(0.5)
newtonLineSearch:set_accept_best(true)
newtonLineSearch:set_check_all(false)
newtonLineSearch:set_suff_descent_factor(0.25)
newtonLineSearch:set_verbose(true)



newtonSolver = NewtonSolver()
newtonSolver:set_linear_solver(solver)
newtonSolver:set_convergence_check(newtonConvCheck)
newtonSolver:set_line_search(newtonLineSearch)
newtonSolver:set_reassemble_J_freq(0)

------------------------------------------------------------------------------------------
-- Apply the solver
------------------------------------------------------------------------------------------




	Interpolate("StartValueX", u, "u")
	Interpolate("StartValueY", u, "v")

	Interpolate("StartValueP", u, "p")
	Interpolate("InitialValue_FractionVolume", u, "c")
	
	
	--LoadVector(u,"currentSolution.vec")
	--print("Reading   ############################")

------------------------------------------------------------------------------------------
-- Prepare the initial guess for the pressure
------------------------------------------------------------------------------------------
if StatBool then
	-- grid function for the solution

	-- Fix the mass fraction and solve the linear problem for the pressure

	fixer = DirichletBoundary()
	domainDisc:add(fixer)
	fixer:invert_subset_selection()
	fixer:add("c", "")

	newtonSolver:init(AssembledOperator(domainDisc))
	newtonSolver:prepare(u)

	-- apply the solver for the stationary pressure problem
	if not newtonSolver:apply(u) then
		print("===> THE PREPARATION PHASE FAILED! <===")
		exit()
	end

	domainDisc:remove (fixer)
end


------------------------------------------------------------------------------------------
-- Apply the solver
------------------------------------------------------------------------------------------
-- start
time = 0
step = 0


	-- compute initial vorticity
	--vorticity(vort,u)

	-- write start solution
	print("Writing start values")
	out = VTKOutput()
	out:clear_selection()
	out:select_all(false)
	out:select_nodal ("u,v", "velocity")
	out:select_nodal ("u", "u")
	out:select_nodal ("v", "v")
	out:select_nodal ("p", "p")
	out:select_nodal ("c", "c")
	out:select(Density, "Rho")
	out:select(Visc, "Mu")


	out:print_subsets(vtk_file_name, u,allSubsets,0,0)
	
newtonSolver:init(op)


if newtonSolver:prepare(u) == false then
	print ("Newton solver prepare failed.") exit()
end
    
-- create new grid function for old value
uOld = u:clone()

tBefore = os.clock()

-- store grid function in vector of  old solutions
solTimeSeries = SolutionTimeSeries()
solTimeSeries:push(uOld, time)


s=1
N_steps=0
for step = 1, numTimeSteps do
	print("++++++ TIMESTEP " .. step .. " BEGIN ++++++")
	
	NewtonSolution=false
	while NewtonSolution==false  do

		-- choose time step
		do_dt = dt
		
		-- setup time Disc for old solutions and timestep
		timeDisc:prepare_step(solTimeSeries, do_dt)
	
		-- prepare newton solver
		if newtonSolver:prepare(u) == false then 
			print ("Newton solver failed at step "..step.."."); exit(); 
		end 
	
		-- apply newton solver
		
		NewtonSolution=false

		if newtonSolver:apply(u)  == false then 
				print ("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<            Reducing TImestep in time, step  " .. step .. " with time step " .. do_dt)
			dt = math.max(dt*red_factor,0.99999*DT_min)
			print("DT=" .. dt .. "");
			print("Time=" .. time .. "");
			if dt < DT_min  or modifyDT== false then 
				print ("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx            Time step below minimum. Aborting. Failed at point with time step" .. do_dt .. "."); exit(); 
	

			else
				do_dt = dt
				VecScaleAssign(u, 1.0, solTimeSeries:latest())
			
			end
			N_steps=0;
			s=0;
		else
			NewtonSolution=true
				
		end 
			
	end

	-- update new time
	time= timeDisc:future_time()
		
	-- get oldest solution
	oldestSol = solTimeSeries:oldest()

	-- copy values into oldest solution (we reuse the memory here)
	VecAssign(oldestSol, u)
	
	-- push oldest solutions with new values to front, oldest sol pointer is poped from end
	solTimeSeries:push_discard_oldest(oldestSol, time)
	
	
	-- compute CFL number
	 
	CFL=cflNumber(u,do_dt)
	print("DT=" .. dt .. "");
	print("Time=" .. time .. "");
	N_steps=N_steps+1
	if(N_steps>=UpdateDt) then
		if modifyDT then dt=math.min(incr_factor*dt,DT_max)  end
		print ("-------------------------------------------------------------------------------------------------Time step increased at Step " .. step .. ", dt =  " .. dt .. ". ")
		N_steps=s*UpdateDt
		newtonConvCheck:set_maximum_steps(max_newton_steps_2)
		newtonSolver:set_line_search(newtonLineSearch)
	end
		
	-- compute kinetic energy
	--ke=kineticEnergy(u)
	--writeNumbers("kineticEnergy.m",step+1,time,ke)
	
	if step % outputFactor == 0 then
	
		out = VTKOutput()
		out:clear_selection()
		out:select_all(false)
		out:select_nodal("u,v", "velocity")
		out:select_nodal("u", "u")
		out:select_nodal("v", "v")
		out:select_nodal("p", "p")
		out:select_nodal("c", "c")
		out:select(Density, "Rho")
		out:select(Visc, "Mu")
		

		out:print_subsets(vtk_file_name, u,allSubsets,step,time)
		print(" ")
	end
	print("++++++ TIMESTEP " .. step .. "  END ++++++")
end

tAfter = os.clock()
newtonSolver:print_average_convergence()
print("Computation took " .. tAfter-tBefore .. " seconds.")

print("done.")
