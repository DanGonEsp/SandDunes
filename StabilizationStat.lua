------------------------------------------------------------------------------------------
-- Navier-Stokes equation, 2d
-- Discretization: Vertex-centered, stabilized
------------------------------------------------------------------------------------------

-- Load utility scripts (e.g. from from ugcore/scripts)
ug_load_script ("ug_util.lua")
ug_load_script ("util/load_balancing_util.lua")

------------------------------------------------------------------------------------------
-- Get command line parameters
------------------------------------------------------------------------------------------

-- Physical parameters
dim 		= util.GetParamNumber("-dim", 2, "dimensionality of the problem")
geometry	= util.GetParam ("-geom", "Dune2D_tri_5")
nu_a 	= util.GetParamNumber("-visc_a", 1.48e-05, "kinematic viscosity")
nu_s 	= util.GetParamNumber("-visc_s", 1.48e-05, "kinematic viscosity")
rho_a 	= util.GetParamNumber("-rho_a", 1.2, "Air Density")
rho_s 	= util.GetParamNumber("-rho_s", 2000, "Sand Density")
dp 	= util.GetParamNumber("-diameter", 0.5e-03, "Particle Diameter")
inflow		= util.GetParamNumber("-inflow", 5, "max. inflow velocity")
c_init		= util.GetParamNumber("-initial concentration", 1.0, "max volume fraction")


alpha_max		= util.GetParamNumber("-max_concentration", 0.6, "max volume fraction")
packing_factor		= util.GetParamNumber("-packing_factor", 0.599, "max volume fraction")
alpha_min		= util.GetParamNumber("-min concentration", 0.5, "max volume fraction")
viscosity_model		= util.GetParamNumber("-granular_model", 2, "Options:  0 Constant, 1 Proportional, 2 Einstein model, 3 Rheology(I), 4 Rheology(I) + Einstein model")
density_model  = util.GetParam("-density_model", "linear", "constant, linear")
interface_value  = util.GetParamNumber("-interface_value", -0.5, "interface value")
bStokes 	= util.GetParamNumber("-Stokes", true ,"If defined, only Stokes Eq. computed")

granular_pressure = true
jumpPressure = false
-- Numerical parameters
harmonic=false
interface_harmonic=false

-- Numerical parameters of the discretization
numRefs 	= util.GetParamNumber("-numRefs", 2, "number of grid refinements")
numPreRefs 	= util.GetParamNumber("-numPreRefs", 2, "number of prerefinements (parallel)")
bNoLaplace 	= util.GetParamNumber("-noLaplace", false,"If defined, only laplace term used")
bExactJac 	= util.GetParamNumber("-exactJac", false,"If defined, exact jacobian used")
bPecletBlend= util.GetParamNumber("-PecletBlend", false,"If defined, Peclet Blend used")
upwind      = util.GetParam("-upwind", "full", "Upwind type full or lps")
bPac        = util.GetParamNumber("-pac", false,"If defined, pac upwind used")
diffLength  = util.GetParam("-difflength", "cor", "fivepoint, raw, corDiffusion length type")

if (jumpPressure) then
	file_name ="PressureStationary"
	stab        = util.GetParam("-stab", "fields_2", "Stabilization type (fields or flow viscosity or karimian)")
else
	file_name ="NoPressureStationary2"
	stab        = util.GetParam("-stab", "fields_2", "Stabilization type (fields or flow viscosity or karimian)")
end

-- Parameters of the solver
ilu_beta	= util.GetParamNumber("-iluBeta", -0.1, "choose a negative value depending on the convection rate")
linIter = util.GetParamNumber("-linIter", 2000, "Max number of linear iterations")
-- Grid file name
gridName = geometry .. ".ugx"

-- Subsets used in the problem
--allSubsets = "Inner,Inlet, Outlet,TopWall, BottomWall, RightWall, Dune,Dune1"
allSubsets = "Inner,Left, Right,Top, Bottom"
--allSubsets = "Inner," .. walls .. ", Inlet, Outlet"

print (" Geometry: " .. geometry .. " (file " .. gridName .. "), dim = " .. dim)
print (" Physical parameter:")
print ("	visc		= " .. nu_a)
print ("	inflow		= " .. inflow)

print (" Numerical parameter:")
print ("	numRefs		= " .. numRefs)
print ("	numPreRefs	= " .. numPreRefs)
print ("	Stokes		= " .. tostring (bStokes))
print ("	noLaplace	= " .. tostring (bNoLaplace))
print ("	exactJac	= " .. tostring (bExactJac))
print ("	PecletBlend = " .. tostring (bPecletBlend))
print ("	upwind		= " .. upwind)
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
			return inflow *0.0
		else
			return inflow *0.0
		end
	else
		return inflow*0.0
	end
end]]
function StartValueX(x,y) 
	hh=14.1856
	nn=2.5
	cc=1.0--math.pow(y/hh,5)
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
	k=0.5/dx
	
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
		value= VolumeFraction2(x,y) 
	--else
		--value = BoundaryVolumeFraction(x,y)
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

---------------------------------------------------------------------- Density


Density = GranularDensityLinker(); 
Density:set_fluid_density(rho_a)
Density:set_particle_density(rho_s*packing_factor)
Density:set_harmonic_value(harmonic)
Density:set_surface_value(interface_harmonic)
Density:set_model(density_model)
Density:set_interface_volume_fraction(interface_value)
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
Visc:set_harmonic_value(harmonic)
Visc:set_surface_value(interface_harmonic)
Visc:set_interface_volume_fraction(interface_value)
Visc:set_limit(1e+02)


PjumpShape= JumpShapeLinker()
PjumpShape:set_interface_volume_fraction(interface_value)


normal = InterfaceNormalLinker()
normal:set_interface_volume_fraction(interface_value)


Ps = ParticlePressureLinker(); 
Ps:set_particle_diameter(dp)
Ps:set_particle_density(rho_s)
Ps:set_fluid_Visc(nu_a*rho_a)
Ps:set_alpha_max(alpha_max)
Ps:set_alpha_min(alpha_min)
Ps:set_packing_factor(packing_factor)

------------------------------------------------------------------------------------------
-- Compose the discretization
------------------------------------------------------------------------------------------

-- inner space

NavierStokesDisc = NavierStokesFV1 (fct_cmp_tbl, {"Inner"})
NavierStokesDisc:set_exact_jacobian (bExactJac)
NavierStokesDisc:set_stokes (bStokes)
NavierStokesDisc:set_laplace ( bNoLaplace)
NavierStokesDisc:set_kinematic_viscosity (Visc)
NavierStokesDisc:set_density(Density)


NavierStokesDisc:set_upwind (upwind)
NavierStokesDisc:set_peclet_blend (bPecletBlend)
NavierStokesDisc:set_stabilization (stab, diffLength)
NavierStokesDisc:set_pac_upwind (bPac)

if (jumpPressure) then
	VolFraction = GridFunctionNumberData(u, "c");
	NavierStokesDisc:set_jump_shape(PjumpShape,diffLength)
	NavierStokesDisc:set_interface_normal(normal)
	NavierStokesDisc:set_vol_fraction(VolFraction)
	NavierStokesDisc:set_interface_value(interface_value)
end


InletDisc = NavierStokesInflow (NavierStokesDisc)
--InletDisc:add ("inflowVel3d", "Inlet")
InletDisc:add ("inflowVel2d", "Left,Top")


SymDiscTop=NavierStokesSymBCFV1(NavierStokesDisc)
SymDiscTop:add("Top")
--SlipDiscTop=NavierStokesWSBCFV1(NavierStokesDisc)
--SlipDiscTop:add("Top")
--SlipDiscTop:set_sliding_factor(0.0000001)
--SlipDiscTop:set_sliding_limit(0.0)

--InletDisc:add ("inflowVel3d", "TopWall")
--InletDisc:add ("inflowVel3d", "RightWall")
--InletDisc:add ("inflowVel3d", "LeftWall")
--FixPressureDisc = DirichletBoundary()
--FixPressureDisc:add(0, "p", "Outlet")

-- boundary condition at the outlet
OutletDisc = NavierStokesNoNormalStressOutflow (NavierStokesDisc)
OutletDisc:add ("Right")



-- boundary condition at the impermeable walls
WallDisc = NavierStokesWall (NavierStokesDisc)
WallDisc:add ("Bottom")




TransportEq = ConvectionDiffusion("c", "Inner", "fv1")
--TransportEq:set_velocity(NavierStokesDisc:velocity_ip())
--TransportEq:set_velocity(NavierStokesDisc:velocity())
TransportEq:set_diffusion(0)
TransportEq:set_upwind(UpwindFV1(upwind)) --upwind type for the transport equation: "no", "full" or "partial"


-- create dirichlet boundary for concentration
dirichletBND = DirichletBoundary()
dirichletBND:add("BoundaryVolumeFraction", "c", "Bottom")
--dirichletBND:add(0.001, "c", "Bottom")

OutflowBND = ConvectionDiffusionOutflowFV1(TransportEq)
OutflowBND:add( "Right, Top, Left, Bottom")	

NeumannBND = NeumannBoundaryFV1("c")
--NeumannBND:add( BottomFlux,"Bottom,Top,Left","Inner")	
--NeumannBND:add( Flux,"Right","Inner")

fixer = DirichletBoundary()	
fixer:invert_subset_selection()
fixer:add("c", "")



print("Transport Equation created.")

---------------------------------------------------------------------------------------
-- Viscosity, Density and Gravitation Input
---------------------------------------------------------------------------------------


Visc:set_volume_fraction(TransportEq:value())
Visc:set_velocity_gradient(NavierStokesDisc:velocity_grad())
if (granular_pressure) then
	Visc:set_particle_pressure(Ps)
	Ps:set_volume_fraction(TransportEq:value())
	Ps:set_velocity_gradient(NavierStokesDisc:velocity_grad())
else
	Visc:set_particle_pressure(NavierStokesDisc:pressure())
end

Density:set_volume_fraction(TransportEq:value())


normal:set_volume_fraction(TransportEq:value())
normal:set_volume_grad(TransportEq:gradient())

PjumpShape:set_volume_fraction(TransportEq:value())




-- the global discretization
domainDisc = DomainDiscretization (approxSpace)
domainDisc:add (NavierStokesDisc)
domainDisc:add (InletDisc)
domainDisc:add (OutletDisc)
domainDisc:add (WallDisc)
--domainDisc:add (SymDiscTop)
--domainDisc:add(FixPressureDisc)
--domainDisc:add (SlipDiscTop)



domainDisc:add(TransportEq)
domainDisc:add(fixer)



-- create the assembled operator for the solver
assembledOp = AssembledOperator (domainDisc)

------------------------------------------------------------------------------------------
-- Set up the solver
------------------------------------------------------------------------------------------

-- For debugging only (to write the intermediate data): --
--util.debug_dir = "FLOW_DEBUG"
--util.debug = { vtk = true, conn_viewer = true }
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
			rap = false,
			smoother =
			{
				type = "ilu",
				beta = ilu_beta,
				damping 	= 0.99,
				--sort	= false,
				--sortEps 	= 1.e-50,
				inversionEps 	= 1.e-24,
				consistentInterfaces   = true
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
			reduction	= 0.5e-8,
			verbose		= true
		}
	},
	lineSearch =
	{
		type			= "standard",
		maxSteps		=7,
		lambdaStart		= 2.0,
		lambdaReduce	= 0.5,
		acceptBest 		= true,
		checkAll		= false
	},
	convCheck =
	{
		type		= "standard",
		iterations	= 200,
		absolute	= 1e-10,
		reduction	= 1e-12,
		verbose		= true
	}
}

solver = util.solver.CreateSolver(solverDesc)

------------------------------------------------------------------------------------------
-- Apply the solver
------------------------------------------------------------------------------------------


Interpolate(StartValueX, u, "u")
Interpolate(StartValueY, u, "v")

Interpolate(StartValueP, u, "p")
Interpolate("InitialValue_FractionVolume", u, "c")

-- Order the DoFs:
OrderLex (approxSpace,  "x")
--OrderCuthillMcKee(approxSpace,true)

out = VTKOutput()
out:clear_selection ()
out:select_nodal (vel_cmp_tbl, "vel")
out:select_nodal ("u", "vel_u")
out:select_nodal ("v", "vel_v")
if dim == 3 then
	out:select_nodal ("w", "vel_w")
end
out:select_nodal ("p", "p")
out:select_nodal ("c", "c")
out:select(Density, "Rho")
out:select(Visc, "Mu")
out:print_subsets ("Test", u, "Inner,Left,Right,Top,Bottom")


-- initialize the solver
solver:init (assembledOp)
solver:prepare (u)

-- apply the solver
if not solver:apply (u) then
	print ("The solver failed.")
	exit ()
end
--vorticity(vort,u)
------------------------------------------------------------------------------------------
-- Plot the results
------------------------------------------------------------------------------------------

out = VTKOutput()
out:clear_selection ()
out:select_nodal (vel_cmp_tbl, "vel")
out:select_nodal ("u", "vel_u")
out:select_nodal ("v", "vel_v")
if dim == 3 then
	out:select_nodal ("w", "vel_w")
end
out:select_nodal ("p", "p")
out:select_nodal ("c", "c")
out:select(Density, "Rho")
out:select(Visc, "Mu")



vtk_file_name = file_name .. "-lev" .. numRefs
if bStokes then
	vtk_file_name = vtk_file_name .. "-Stokes"
end
print ("Output to file " .. vtk_file_name .. ".vtu")
out:print_subsets (vtk_file_name, u, "Inner,Left,Right,Top,Bottom")


--"Inner,UpperWall,LowerWall,FrontWall,BackWall,CylinderWall,Inlet,Outlet")

-- End of File
