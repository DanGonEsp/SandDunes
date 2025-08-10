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

-- Geometry parameters
geom_name = "quad"
if geom_name == "tri" then
	geometry	= util.GetParam ("-geom", "Dune2D_tri_5")
	file_name = "Tri"
else
	geometry	= util.GetParam ("-geom", "Dune2D_quads")
	file_name = "Quads"
end
gridName = geometry .. ".ugx"
-- Subsets used in the problem
allSubsets = "Inner,Left, Right,Top, Bottom"


-- Physical parameters
dim 		= util.GetParamNumber("-dim", 2, "dimensionality of the problem")
nu_a 	= util.GetParamNumber("-visc_a", 1.48e-05, "kinematic viscosity")
rho_a 	= util.GetParamNumber("-rho_a", 1.2, "Air Density")
rho_s 	= util.GetParamNumber("-rho_s", 2500, "Sand Density")
dp 	= util.GetParamNumber("-diameter", 1e-03, "Particle Diameter")
nu_s 	= util.GetParamNumber("-visc_s", 1.48e-05, "kinematic viscosity")
inflow		= util.GetParamNumber("-inflow", 10, "max. inflow velocity")
c_init		= util.GetParamNumber("-initial concentration", 1.0, "max volume fraction")


alpha_max		= util.GetParamNumber("-max_concentration", 0.635, "max volume fraction")
packing_factor		= util.GetParamNumber("-packing_factor", 0.625, "max volume fraction")
alpha_min		= util.GetParamNumber("-min concentration", 0.57, "max volume fraction")
viscosity_model		= util.GetParamNumber("-granular_model", 3, "Options:  0 Constant, 1 Proportional, 2 Einstein model, 3 Rheology(I) + Einstein model")
density_model  = util.GetParam("-density_model", "linear", "constant, linear")
interface_value  = util.GetParamNumber("-interface_value",  alpha_min*c_init/packing_factor, "interface value")
--interface_value  = util.GetParamNumber("-interface_value",  0.3, "interface value")
FR = 0.05
B_phi = 1
deltaGamma = 1e-04;
Visc_limit = 1e6

deltaPs = 1.48e-04;
deltaI = 1e-03;
FricMu_1=0.38
FricMu_2=0.64
I_0 = 0.279




jumpPressure = false
boolSource = true
StatBool = true


bodyForceSource = false
boolGradientPsSource = false


if jumpPressure then

	file_name =file_name .. "Pressure"
	value_beta =-0.005         --value_beta = -0.005
	jac = 0.0
	damping_mg = 1.0        --damping_mg = 0.15
else
	file_name =file_name .. "NoPressure"
	value_beta = -0.5
	jac = false
	damping_mg = 1.0
end
stab        = util.GetParam("-stab", "fields_2", "Stabilization type (fields or flow viscosity or karimian)")

-- Numerical parameters of the discretization
numRefs 	= util.GetParamNumber("-numRefs", 2, "number of grid refinements")
numPreRefs 	= util.GetParamNumber("-numPreRefs", 1, "number of prerefinements (parallel)")

bStokes 	= util.GetParamNumber("-Stokes", false ,"If defined, only Stokes Eq. computed")
bNoLaplace 	= util.GetParamNumber("-noLaplace", false,"If defined, only laplace term used")
bExactJac 	= util.GetParamNumber("-exactJac", jac,"If defined, exact jacobian used")
bPecletBlend= util.GetParamNumber("-PecletBlend", false,"If defined, Peclet Blend used")
upwind      = util.GetParam("-upwind", "full", "Upwind type full or lps")
bPac        = util.GetParamNumber("-pac", false,"If defined, pac upwind used")
diffLength  = util.GetParam("-difflength", "cor", "fivepoint, raw, corDiffusion length type")


dt = util.GetParamNumber("-dt",1.0)
numTimeSteps =  util.GetParamNumber("-numTimeSteps", 20	)
EndTime = util.GetParamNumber("-EndTime", 6.2, "EndTime")
boolEndTime 	= util.GetParamNumber("-boolEndTime", true)
outputFactor     = util.GetParam("-output", 1, "output every ... steps")
max_newton_steps=util.GetParamNumber("-numNewtonSteps", 200)
max_linear_steps=util.GetParamNumber("-numLinearIter", 200)
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



-- Parameters of the solver
linIter = util.GetParamNumber("-linIter", 2000, "Max number of linear iterations")










vtk_file_name = file_name .. "-lev" .. numRefs
if boolSource then
	if bodyForceSource then
		vtk_file_name = vtk_file_name .. "-BodyForce"
	else
		vtk_file_name = vtk_file_name .. "-GradientForce"
	end
	if boolGradientPsSource then
		vtk_file_name = vtk_file_name .. "-PartPressGradient"
	else
		vtk_file_name = vtk_file_name .. "-NoPartPressGradient"
	end
else
	vtk_file_name = vtk_file_name .. "-NoForce"
end
if bStokes then
	vtk_file_name = vtk_file_name .. "-Stokes"
end




print (" Geometry: " .. geometry .. " (file " .. gridName .. "), dim = " .. dim)
print (" Physical parameter:")
print ("	inflow		= " .. inflow)
print ("	Stokes		= " .. tostring (bStokes))
print ("	Pressure Jump	= " .. tostring (jumpPressure))
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
if dim == 3 then
	approxSpace:add_fct("w", "Lagrange",1,allSubsets)
end
approxSpace:add_fct("p", "Lagrange",1,allSubsets)
approxSpace:add_fct("c", "Lagrange",1,allSubsets)



approxSpace:init_levels()
approxSpace:init_top_surface()
approxSpace:print_statistic()

util.solver.defaults.approxSpace = approxSpace


-- grid function for the solution
u = GridFunction (approxSpace)
u:set(0)

--------------------------------
--------------------------------
-- Lua Functions
--------------------------------
--------------------------------
------------------------------------------------------------------------------------------
-- Lua Functions
------------------------------------------------------------------------------------------
H_domain=5
H_0=2
h_0=0.0
sigma=3--0.4
mu_c=20
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
k2=0.05 
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
function StartValueX(x,y) 
	hh=14.1856
	nn=2.5
	cc=math.pow(y/hh,5)
	return inflow*(math.min(1.0, math.pow(y/hh,1/nn))*(1-cc) +(cc)* (2*hh - y) * (y ) / (hh * hh))
end

function StartValueY(x,y) 
	return 0.0*inflow
end
function StartValueZ(x,y) 
	return 0
end

---------------------------------------------------------------------- Initial Pressure
Pstd=0.0	
function StartValueP(x,y) 
	if y>Dune(x) then
		return Pstd *(x+30)
	else 
		return  Pstd *(x+30)
	end
end
function PressureBoundary(x,y) 
	return  -1.2*y*9.81
end

---------------------------------------------------------------------- Initial VolumeFraction 



--[[function VolumeFraction(x,y)
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
	
	
end]]
--[[function VolumeFraction(x,y)

	if y>5 then
		return  0.0
	else
		return  c_init
	end
end]]
function VolumeFraction2(x,y)
	dd=dist(x,y)
	ds=0*dx
	kk1=1600
	kk2=1600 
	
	if y>Dune(x,y) then
		if dd<ds then
			dd=-dd
			k=k
			return c_init / (1.0 + math.exp(-kk1*dd))
		else
			return 0.0
		end
	else
		if dd<ds then
			k=k
			return c_init / (1.0 + math.exp(-kk2*dd))
		else
			return c_init 
		end
	end
	
	
	
end
function InitialValue_FractionVolume(x,y)
	--if y>0 then 
		value= VolumeFraction2(x,y) 
	--else
		--value = VolumeFraction(x,y)
	--end
	return value
end
---------------------------------------------------------------------- Boundary Condition  
----------------------------------------------------------- Inlet
local H=1

function inflowVel2d(x, y, t)

	return StartValueX(x,y),StartValueY(x,y)
end

----------------------------------------------------------- Bottom

function BoundaryVolumeFraction(x,y)
	a= 10
	b = 20
	if x>a and x<b then--y>Dune(x,y) then
		return c_init
	else
		return 0.0
	end
end

function BottomFlux(x,y) return 0 end

---------------------------------------------------------------------- VolFraction
ccc = GridFunctionNumberData(u, "c");
VolFraction = ScaleLinker();
VolFraction:set_import_1(ccc)
VolFraction:set_import_2(1.0)

---------------------------------------------------------------------- Density

Density = GranularDensityLinker(); 
Density:set_fluid_density(rho_a)
Density:set_particle_density(rho_s*packing_factor)
Density:set_model(density_model)
Density:set_interface_volume_fraction(interface_value)
---------------------------------------------------------------------- Particle Pressure


Ps = ParticlePressureLinker(); 
Ps:set_particle_diameter(dp)
Ps:set_particle_density(rho_s)
Ps:set_fluid_Visc(nu_a*rho_a)
Ps:set_alpha_max(alpha_max)
Ps:set_alpha_min(alpha_min)
Ps:set_packing_factor(packing_factor)
Ps:set_mix_density(Density)
Ps:set_gravity(-9.81)
Ps:set_FR(FR)
Ps:set_B_phi(B_phi)
Ps:set_deltaGamma(deltaGamma)

---------------------------------------------------------------------- Viscosity
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
Visc:set_interface_volume_fraction(interface_value)
Visc:set_limit(Visc_limit)
Visc:set_particle_pressure(Ps)
Visc:set_deltaPs(deltaPs)
Visc:set_deltaI(deltaI)
Visc:set_FricMu_1(FricMu_1)
Visc:set_FricMu_2(FricMu_2)
Visc:set_I_0(I_0)
Visc:set_deltaGamma(deltaGamma)



PjumpShape= JumpShapeLinker()
PjumpShape:set_interface_volume_fraction(interface_value)


normal = InterfaceNormalLinker()
normal:set_interface_volume_fraction(interface_value)

GradientPsSource = ParticlePressureGradientLinker()
GradientPsSource:set_particle_pressure(Ps)
GradientPsSource:set_interface_volume_fraction(interface_value)

Source = GranularSourceLinker()
Source:set_particle_density(rho_s)
Source:set_fluid_density(rho_a)
Source:set_mix_density(Density)
Source:set_gravity(-9.81)
Source:set_packing_factor(packing_factor)
Source:set_rel_vel(0.0)
Source:set_interface_volume_fraction(interface_value)
Source:set_bool_body_force(bodyForceSource)
Source:set_bool_particle_pressure_force(boolGradientPsSource)
Source:set_ps_grad(GradientPsSource)

------------------------------------------------------------------------------------------
-- Compose the discretization
------------------------------------------------------------------------------------------

-- inner space

NavierStokesDisc = NavierStokesFV1 (fct_cmp_tbl, {"Inner"})
NavierStokesDisc:set_exact_jacobian (bExactJac)
NavierStokesDisc:set_stokes (bStokes)
NavierStokesDisc:set_laplace ( bNoLaplace)
NavierStokesDisc:set_upwind (upwind)
NavierStokesDisc:set_peclet_blend (bPecletBlend)
NavierStokesDisc:set_stabilization (stab, diffLength)
NavierStokesDisc:set_pac_upwind (bPac)

NavierStokesDisc:set_kinematic_viscosity (Visc)
NavierStokesDisc:set_density(Density)
NavierStokesDisc:set_density_ref(0.0)
if (boolSource) then
	NavierStokesDisc:set_source(Source)
end
NavierStokesDisc:set_interface_value(interface_value)


InletDisc = NavierStokesInflow (NavierStokesDisc)
--InletDisc:add ("inflowVel3d", "Inlet,Top,Bottom")
InletDisc:add ("inflowVel2d","Left,Top")

-- boundary condition at the outlet
OutletDisc = NavierStokesNoNormalStressOutflow (NavierStokesDisc)
OutletDisc:add ("Right")

-- boundary condition at the impermeable walls
WallDisc = NavierStokesWall (NavierStokesDisc)
WallDisc:add ("Bottom")

--PressureOutlet = DirichletBoundary()
--PressureOutlet:add("PressureBoundary", "p", "Bottom")




Vel = VelocityBCLinker()
--Stress = NavierStokesInflowStressFV1(NavierStokesDisc)
--Stress:add("Left,Top")
--Stress:set_velocity(Vel)





TransportEq = ConvectionDiffusion("c", "Inner", "fv1")
TransportEq:set_velocity(NavierStokesDisc:velocity_ip())
--TransportEq:set_velocity(NavierStokesDisc:velocity())
TransportEq:set_diffusion(0)
TransportEq:set_upwind(UpwindFV1(upwind)) --upwind type for the transport equation: "no", "full" or "partial"
TransportEq:set_interface_value(interface_value)


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

Density:set_volume_fraction(TransportEq:value())

Visc:set_volume_fraction(TransportEq:value())
Visc:set_velocity_gradient(NavierStokesDisc:velocity_grad())


PjumpShape:set_volume_fraction(TransportEq:value())

Ps:set_volume_fraction(TransportEq:value())
Ps:set_velocity_gradient(NavierStokesDisc:velocity_grad())


if (jumpPressure) then
	NavierStokesDisc:set_jump_shape(PjumpShape,diffLength)
	NavierStokesDisc:set_interface_normal(normal)
	NavierStokesDisc:set_vol_fraction(VolFraction)
end

normal:set_volume_fraction(TransportEq:value())
normal:set_volume_grad(TransportEq:gradient())

Source:set_volume_fraction(TransportEq:value())
Source:set_volume_grad(TransportEq:gradient())

Vel:set_volume_fraction(TransportEq:value())
Vel:set_volume_grad(TransportEq:gradient())
Vel:set_particle_density(rho_s)
Vel:set_fluid_density(rho_a)
Vel:set_mix_density(Density)
Vel:set_gravity(-9.81)
Vel:set_packing_factor(packing_factor)
Vel:set_rel_vel(0.0)
Vel:set_interface_volume_fraction(interface_value)



-- the global discretization
domainDisc = DomainDiscretization (approxSpace)
domainDisc:add (NavierStokesDisc)
domainDisc:add (InletDisc)
domainDisc:add (OutletDisc)
domainDisc:add (WallDisc)
--domainDisc:add (PressureOutlet)



domainDisc:add(TransportEq)
domainDisc:add(OutflowBND)




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
--util.debug = { vtk = true, conn_viewer = false }
--util.CreateGridFuncDebugWriter (approxSpace)
-- --

solverDesc =
{
	type = "newton",
	debug = false, -- for the debug output from the Newton's method
	linSolver =
	{
		type = "bicgstab",
		precond =
		{
			type = "gmg",
			rap = true,
			smoother =
			{
				type = "ilu",
				beta = value_beta,
				damping 	= damping_mg,
				--sort	= false,
				--sortEps 	= 1.e-50,
				inversionEps 	= 1.e-8,
				--consistentInterfaces   = true
				--overlap 		= false,
				--ordering 		= nil
			},
			preSmooth = 1,
			postSmooth = 1,
			baseSolver = "lu",
			baseLevel = numPreRefs
		},
		convCheck =
		{
			type		= "standard",
			iterations	= linIter,
			absolute	= 1e-12,
			reduction	= 1e-2,
			verbose		= true
		}
	},
	lineSearch =
	{
		type			= "standard",
		maxSteps		=5,
		lambdaStart		= 2.0,
		lambdaReduce	= 0.5,
		acceptBest 		= true,
		checkAll		= false
	},
	convCheck =
	{
		type		= "standard",
		iterations	= max_newton_steps,
		absolute	= 1e-5,
		reduction	= 1e-12,
		verbose		= true
	}
}

solver = util.solver.CreateSolver(solverDesc)

------------------------------------------------------------------------------------------
-- Apply the solver
------------------------------------------------------------------------------------------





	
	

	
------------------------------------------------------------------------------------------
-- Prepare the initial guess for the pressure
------------------------------------------------------------------------------------------
if StatBool then

	Interpolate(StartValueX, u, "u")
	Interpolate("StartValueY", u, "v")

	Interpolate("StartValueP", u, "p")
	Interpolate("InitialValue_FractionVolume", u, "c")

	-- grid function for the solution

	-- Fix the mass fraction and solve the linear problem for the pressure

	fixer = DirichletBoundary()
	domainDisc:add(fixer)
	fixer:invert_subset_selection()
	fixer:add("c", "")

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
	domainDisc:remove (fixer)
	
	print("++++++++++++++++++++++++ INITIAL CONDITIONS  (STEADY STATE DONE) ++++++++++++++++++++++++")
else

	Interpolate(0.0*inflow, u, "u")
	Interpolate("StartValueY", u, "v")

	Interpolate("StartValueP", u, "p")
	Interpolate("InitialValue_FractionVolume", u, "c")
	
end


------------------------------------------------------------------------------------------
-- Apply the solver
------------------------------------------------------------------------------------------
-- start
time = 0
step = 0

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
out:select(Ps, "Ps")
out:print_subsets(vtk_file_name, u,allSubsets,0,0)
print ("Output to file " .. vtk_file_name .. ".vtu  in time t = 0")
solver:init(op)



if solver:prepare(u) == false then
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
		if solver:prepare(u) == false then 
			print ("Newton solver failed at step "..step.."."); exit(); 
		end 
	
		-- apply newton solver
		
		NewtonSolution=false

		if solver:apply(u)  == false then 
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
		if modifyDT then 
			dt=math.min(incr_factor*dt,DT_max)  
			print ("-------------------------------------------------------------------------------------------------Time step increased at Step " .. step .. ", dt =  " .. dt .. ". ")
		end
		
		N_steps=s*UpdateDt
	end
		
	
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
		out:select(Ps, "Ps")
		

		out:print_subsets(vtk_file_name, u,allSubsets,step,time)
		print(" ")
	end
	print("++++++ TIMESTEP " .. step .. "  END ++++++")
end

tAfter = os.clock()
--solver:print_average_convergence()
print("-			-")
print("-------------------------------------------------------------------------------")
print("Steady state Computation took " .. tAfter_s-tBefore_s .. " seconds.")
print("Temporal Computation took " .. tAfter-tBefore .. " seconds.")
print("Total Computation took " .. tAfter_s-tBefore_s+tAfter-tBefore .. " seconds.")
print("-------------------------------------------------------------------------------")
print("")
print("")
print ("Output to file " .. vtk_file_name .. ".vtu")
print("done.")

