------------------------------------------------------------------------------------------
-- Navier-Stokes equation, 3d
-- Discretization: Vertex-centered, stabilized
------------------------------------------------------------------------------------------

-- Load utility scripts (e.g. from from ugcore/scripts)
ug_load_script ("ug_util.lua")
ug_load_script ("util/load_balancing_util.lua")

ug_load_script("util/domain_disc_util.lua")
ug_load_script("util/conv_rates_kinetic.lua")

RequiredPlugins({"Limex", "NavierStokes"})

------------------------------------------------------------------------------------------
-- Split communicator
------------------------------------------------------------------------------------------
local numProc         = util.GetParamNumber("-numProc", 1, "Number of temporal processes")

SpaceTimeComm = SpaceTimeCommunicator()
SpaceTimeComm:split(numProc)

local rank = ProcRank()
local rank_t=SpaceTimeComm:get_temporal_rank()

local TemporalSize = SpaceTimeComm:get_temporal_size()
local SpaceSize = SpaceTimeComm:get_spatial_size()

print("TemporalSize = " ..TemporalSize)
print("SpaceSize = " ..SpaceSize)

------------------------------------------------------------------------------------------
-- Input parameter table
------------------------------------------------------------------------------------------
local csvfile = require "simplecsv"
local InValues, num_rows, num_cols = csvfile.read('./Table_in.csv') -- read file csv1.txt to matrix m
if( TemporalSize > num_rows-1) then print ("TemporalSize larger than rows in input parametrs."); exit(); end


inflow = InValues[rank_t+2][1]
H_0 = InValues[rank_t+2][2]
W0 = InValues[rank_t+2][3]

print("Inflow = " ..inflow.."m/s")
print("Heigh = " ..H_0.. "m.")
print("Width = " ..W0.. "m.")

local fixedNum = string.format("%04d", rank_t)


------------------------------------------------------------------------------------------
-- parameters
------------------------------------------------------------------------------------------
params =
{
			-- Numerical parameters of the discretization
	dim      = util.GetParamNumber("-dim", 2, "dimensionality of the problem"),
	file_name = util.GetParam("-file_name", "Test2") .. fixedNum,
	folder_name = util.GetParam("-folder_name", "Test2"),
	elem_type = util.GetParam("-elem_type", "tri", "tri, quad"),
	numRefs     = util.GetParamNumber("-numRefs", 1, "number of grid refinements"),
	numPreRefs     = util.GetParamNumber("-numPreRefs", 0, "number of prerefinements (parallel)"),
	startTime  = util.GetParamNumber("-start", 0.0, "start time"),
	endTime    = util.GetParamNumber("-end", 1000.0, "end time"),
	numTimeSteps    = util.GetParamNumber("-numTimeSteps", 2, "time steps"),
	DTmin= util.GetParamNumber("-DTmin", 1e-05, "min  DT"),
	outputFactor     = util.GetParam("-output", 1, "output every ... steps"),
	
	timeMethod = util.GetParam("-timeMethod","euler","euler limex"),
	modifyDT     = util.GetParamBool("-modifyDT", false),
	incr_factor     = util.GetParamNumber("-incr_factor", 1.3),
	red_factor_fail     = util.GetParamNumber("-red_factor_fail", 0.7),
	red_factor_success     = util.GetParamNumber("-red_factor_success", 0.97),
	optimal_newton_steps = util.GetParamNumber("-optimal_newton_steps", 20),
	maxConvRate = util.GetParamNumber("-maxConvRate", 0.9),
	minConvRate = util.GetParamNumber("-minConvRate", 0.1),
	boolDebug = util.GetParamBool("-boolDebug", false),
	boolDebugStep = util.GetParamBool("-boolDebugStep", false),
	
	tol     = util.GetParamNumber("-limex-tol", 1e-2, "time step size"),
	nstages = util.GetParamNumber("-limex-nstages", 2, "limex stages (2 default)"),
	limex_partial_mask = util.GetParamNumber("-limex-partial", 0, "limex partial (0 or 3)"),
	limex_debug_level = util.GetParamNumber("-limex-debug-level", 5, "limex debug level (integer)"),
	
	max_newton_steps_steady_state=util.GetParamNumber("-numNewtonSteps", 100),
	max_newton_steps_transcient=util.GetParamNumber("-max_newton_steps_transcient", 50),
	AbsDefect = util.GetParamNumber("-AbsDefect", 1e-05),
	RedDefect = util.GetParamNumber("-RedDefect", 1e-07),

	max_linear_steps=util.GetParamNumber("-max_linear_steps", 200),
	damping_mg = util.GetParamNumber("-damping_mg", 0.9),
	value_beta = util.GetParamNumber("-value_beta", 0.0 ),
	LinAbsDefect = util.GetParamNumber("-LinAbsDefect", 1e-07),
	LinRedDefect = util.GetParamNumber("-LinRedDefect", 1e-02),

	lambdamaxSteps = util.GetParamNumber("-lambdamaxSteps", 5),
	lambdaStart  = util.GetParamNumber("-lambdaStart", 1.0),
	
			-- Physical phenomenon of simulation
	doSteadyState = util.GetParamBool("-doSteadyState", true),
	boolSource = util.GetParamBool("-boolSource", false),
	consistentRho_in_source = util.GetParamBool("-consistentRho_in_source", true),
	boolRelativeVel = util.GetParamBool("-boolRelativeVel", true),
	boolGradientPsSource = util.GetParamBool("-boolGradientPsSource", false),
	boolViscPs = util.GetParamBool("-boolViscPs", true),
	boolAveDiff = util.GetParamBool("-boolAveDiff", true),
	boolSlipDiff = util.GetParamBool("-boolSlipDiff", false),
	boolSlipVel = util.GetParamBool("-boolSlipVel", false),
	
	inflow   = inflow,
	H_0= H_0,
	ReferencePressure  = util.GetParamNumber("-ReferencePressure",  1.7493e2, "interface value"),
	bStokes     = util.GetParamBool("-Stokes", true ,"If defined, only Stokes Eq. computed"),
	bNoLaplace     = util.GetParamNumber("-noLaplace", false,"If defined, only laplace term used"),
	bExactJac     = util.GetParamNumber("-exactJac", 0.0,"If defined, exact jacobian used"),
	bPecletBlend = util.GetParamBool("-PecletBlend", false,"If defined, Peclet Blend used"),
	upwind_m      = util.GetParam("-upwind_m", "full", "Upwind type full or lps"),
	upwind_t      = util.GetParam("-upwind_t", "full", "Upwind type full or lps"),
	upwind_r      = util.GetParam("-upwind_r", "full", "Upwind type full or lps"),
	bPac        = util.GetParamNumber("-pac", false,"If defined, pac upwind used"),
	diffLength  = util.GetParam("-difflength", "raw", "fivepoint, raw, corDiffusion length type"),
	stab        = util.GetParam("-stab", "fields_2", "Stabilization type (fields or flow viscosity or karimian)"),
	turbViscMethod = util.GetParam("-turbulenceModel","no","TurbVismodel type no , dyn or sma"),
	modellconstant = util.GetParamNumber("-c",0.1),

	--Material Properties
	nu_a     = util.GetParamNumber("-visc_a", 1.48e-05, "kinematic viscosity"),
	rho_a     = util.GetParamNumber("-rho_a", 1.2, "Air Density"),
	rho_s     = util.GetParamNumber("-rho_s", 2500, "Sand Density"),
	dp     = util.GetParamNumber("-diameter", 1e-03, "Particle Diameter"),
	nu_s     = util.GetParamNumber("-visc_s", 1.48e-05, "kinematic viscosity"),
	c_init        = util.GetParamNumber("-c_init", 0.6, "max volume fraction"),

	alpha_max        = util.GetParamNumber("-alpha_max", 0.635, "max volume fraction"),
	alpha_min        = util.GetParamNumber("-min alpha_min", 0.57, "max volume fraction"),
	granular_model= util.GetParamNumber("-granular_model", 3, "Opt: 0 Const, 1 Linear, 2 Einstein, 3 Rheology(I) + Einstein"),
	density_model  = util.GetParam("-density_model", "linear", "constant, linear"),
	interface_value  = util.GetParamNumber("-interface_value",  0.57, "interface value"),
	--interface_value  = util.GetParamNumber("-interface_value",  0.3, "interface value"),
	drag_mod = util.GetParamNumber("-drag_model", 2, "Opt: 0 StokesLaw, 1 formula, 2 Schiller-Naumann, 3 Turton and Levenspiel"),
	--Model 0 pow(0.63+4.8/sqrt(RE),2.0);

	FR = 0.05,
	B_phi = 1,
	deltaGamma = 1e-03,
	Visc_limit = 1e15,

	deltaPs = 1.48e-04,
	deltaI = 1e-03,
	FricMu_1=0.38,
	FricMu_2=0.64,
	I_0 = 0.279,
	
}


params.DTmax = (params.endTime - params.startTime) / params.numTimeSteps
dt = util.GetParamNumber("-dt", params.DTmax, "dt")
print("fixedNum = "..fixedNum)
file_name = params.file_name
folder_name = params.folder_name
c_init = params.c_init
dim = params.dim

------------------------------------------------------------------------------------------
-- Get command line parameters
------------------------------------------------------------------------------------------

-- Geometry parameters

if params.elem_type == "tri" then
	gridName	= util.GetParam ("-geom", "Dune2D_tri_double.ugx")
else
	gridName	= util.GetParam ("-geom", "Dune2D_quads_double.ugx")
end
-- Subsets used in the problem
allSubsets = "Inner, Inner2,Left, Right,Top, Bottom"
Inner_total={"Inner","Inner2"}

------------------------------------------------------------------------------------------
-- Folder and files
------------------------------------------------------------------------------------------

if (rank==0) then
	if not DirectoryExists (folder_name) then
		CreateDirectory (folder_name)
	end
	print(folder_name)
end
vtk_file_name = file_name .. "-lev" .. params.numRefs
print(vtk_file_name)
if params.bStokes then
    vtk_file_name = vtk_file_name .. "-Stokes"
end

if params.boolRelativeVel then
    vtk_file_name =vtk_file_name .. "-RelVel"
else
    vtk_file_name =vtk_file_name .. "-NoRelVel"
end


if  params.boolGradientPsSource  or params.boolSource  then

    if params.boolSource then
        if params.consistentRho_in_source then vtk_file_name =vtk_file_name .. "-Consistent" end
        vtk_file_name = vtk_file_name .. "-MG_Force"
    end
    
    if params.boolGradientPsSource then
        vtk_file_name = vtk_file_name .. "-DPs"
    else
        vtk_file_name = vtk_file_name .. "-DPh"
    end
    
else
    vtk_file_name = vtk_file_name .. "-NoForce"
end

if params.boolAveDiff then
    vtk_file_name = vtk_file_name .. "-AveDiff"
else
    vtk_file_name = vtk_file_name .. "-NoAveDiff"
end

if params.boolSlipDiff then
    vtk_file_name = vtk_file_name .. "-SlipDiff"
else if params.boolSlipVel then
		vtk_file_name = vtk_file_name .. "-SlipVel"
	else
		vtk_file_name = vtk_file_name .. "-NoSlip"
	end
end

folder = folder_name .. "/" .. vtk_file_name
if (rank== 0) then
	if not DirectoryExists (folder) then
		CreateDirectory (folder)
	else
		print(" You must erase the content")
	end
end

vtk_file_name = folder .. "/" .. vtk_file_name -- VTK output file name base

------------------------------------------------------------------------------------------
-- Initialize UG4
------------------------------------------------------------------------------------------

InitUG (dim, AlgebraType("CPU", dim+2))
GetLogAssistant():enable_file_output(true, folder .. "/LogFile_"..rank_t.. "_Lev"..params.numRefs)
if rank_t > 0 then GetLogAssistant():enable_terminal_output(false) end

------------------------------------------------------------------------------------------
-- Initialize UG4
------------------------------------------------------------------------------------------
print (" Geometry: " .. gridName ..", dim = " .. params.dim)
print (" Physical parameter:")
print ("	inflow			= " .. params.inflow)
print ("	Stokes			= " .. tostring (params.bStokes))
print ("	Steady state    	= " .. tostring (params.doSteadyState))
print ("	BodyForce       	= " .. tostring (params.boolSource))
print ("	Consisten Gravity	= " .. tostring (params.consistentRho_in_source))
print ("	Relative Vel    	= " .. tostring (params.boolRelativeVel))
print ("	Ps gradient     	= " .. tostring (params.boolGradientPsSource))
print ("	Ps in visc      	= " .. tostring (params.boolViscPs))
print ("	Diffusion       	= " .. tostring (params.boolAveDiff))
print ("	SlipDiff         	= " .. tostring (params.boolSlipDiff))
print ("	SlipVel         	= " .. tostring (params.boolSlipVel))
print (" Numerical parameter:")
print ("	numRefs			= " .. params.numRefs)
print ("	numPreRefs		= " .. params.numPreRefs)
print ("	noLaplace		= " .. tostring (params.bNoLaplace))
print ("	exactJac		= " .. tostring (params.bExactJac))
print ("	PecletBlend		= " .. tostring (params.bPecletBlend))
print ("	upwind_m		= " .. params.upwind_m)
print ("	upwind_t		= " .. params.upwind_t)
print ("	pac			= " .. tostring (params.bPac))
print ("	stab			= " .. params.stab)
print ("	difflength		= " .. params.diffLength)

--------------------------------------------------------------------------------
-- Problem setup.
--------------------------------------------------------------------------------
local myProblem=require("SandDunesConfig")
myProblem:Init(params)

------------------------------------------------------------------------------------------
-- load, refine and distribute the grid
------------------------------------------------------------------------------------------

if dim == 3 then
	fct_cmp_tbl = {"u", "v", "w", "p", "c"}
	vel_cmp_tbl = {"u", "v", "w"}
else
	fct_cmp_tbl = {"u", "v", "p", "c"}
	vel_cmp_tbl = {"u", "v"}
end	

-- Create the domain, load the grid and refine it
dom = util.CreateDomain (gridName, params.numPreRefs)
balancer.RefineAndRebalanceDomain (dom, params.numRefs - params.numPreRefs)
--util.refinement.CreateRegularHierarchy(dom, params.numRefs, true)

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
h_0=0.0
sigma=3
mu_c=20

mu_c1=mu_c-W0/2
mu_c2=mu_c+W0/2
sigma1=sigma
sigma2=0.3*sigma
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

dx=math.pow(1/2,(params.numRefs+2))
---------------------------------------------------------------------- Initial Velocity
--[[function StartValueX(x,y)
	hh=14.1856
	nn=2.5
	cc=math.pow(y/hh,5)
	return inflow*(math.min(1.0, math.pow(y/hh,1/nn))*(1-cc) +(cc)* (2*hh - y) * (y ) / (hh * hh))
end]]
function StartValueX(x,y)
	hh=14.1856
	return inflow* (2*hh - y) * (y ) / (hh * hh)
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




function VolumeFraction(x,y)
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


function StartValueC(x,y)
	return VolumeFraction(x,y)
end
---------------------------------------------------------------------- Boundary Condition
----------------------------------------------------------- Inlet
local H=1

function inflowVel2d(x, y, t)

	return StartValueX(x,y),StartValueY(x,y)
end

function ConstValue(x, y, t)
	if (y < 0.88667) then
		return 1.2e-02
	else return 0.0
	end
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
---------------------------------------------------------------------- Interface Properties

InterfaceValues = Interface()

InterfaceValues:set_particle_diameter(params.dp)
InterfaceValues:set_particle_density(params.rho_s)
InterfaceValues:set_air_density(params.rho_a)
InterfaceValues:set_fluid_Visc(params.nu_a*params.rho_a)
InterfaceValues:set_alpha_max(params.alpha_max)
InterfaceValues:set_packing_factor(params.c_init)
InterfaceValues:set_alpha_min(params.alpha_min)
InterfaceValues:set_FR(params.FR)
InterfaceValues:set_B_phi(params.B_phi)
InterfaceValues:set_deltaGamma(params.deltaGamma)
InterfaceValues:set_limit(params.Visc_limit)
InterfaceValues:set_bool_particle_pressure_force(params.boolGradientPsSource)
InterfaceValues:set_bool_consistent_gravity(false)
InterfaceValues:set_reference_pressure(params.ReferencePressure)
InterfaceValues:set_time_step_factor(dt)
InterfaceValues:set_interface_volume_fraction(params.interface_value)
InterfaceValues:set_drag_model(params.drag_mod)
InterfaceValues:set_bool_initialized(true)

---------------------------------------------------------------------- Density

Density = GranularDensityLinker(); 
Density:set_fluid_density(params.rho_a)
Density:set_particle_density(params.rho_s)
Density:set_model(params.density_model)
Density:set_interface_volume_fraction(params.interface_value)

---------------------------------------------------------------------- Viscosity
Visc = GranularViscosityLinker(); 

Visc:set_granular_model(params.granular_model)
Visc:set_particle_diameter(params.dp)
Visc:set_particle_density(params.rho_s)
Visc:set_particle_kinematicVisc(params.nu_s)
Visc:set_fluid_density(params.rho_a)
Visc:set_fluid_Visc(params.nu_a*params.rho_a)
Visc:set_alpha_max(params.alpha_max)
Visc:set_alpha_min(params.alpha_min)
Visc:set_mix_density(Density)
Visc:set_interface_volume_fraction(params.interface_value)
Visc:set_deltaPs(params.deltaPs)
Visc:set_deltaI(params.deltaI)
Visc:set_FricMu_1(params.FricMu_1)
Visc:set_FricMu_2(params.FricMu_2)
Visc:set_I_0(params.I_0)
Visc:set_deltaGamma(params.deltaGamma)



if params.turbViscMethod=="dyn" then
	viscosityData = FV1DynamicTurbViscData(approxSpace,u)
	
else 
	viscosityData = FV1SmagorinskyTurbViscData(approxSpace,u,params.modellconstant)
 end
viscosityData:set_turbulence_zero_bnd("Left")
viscosityData:set_kinematic_viscosity(Visc)

---------------------------------------------------------------------- Scale

Viscosity = ScaleLinker();
Viscosity:set_import_1(Density)
Viscosity:set_import_2(Visc)

Viscosity2 = ScaleLinker();
Viscosity2:set_import_1(Density)
Viscosity2:set_import_2(viscosityData)
---------------------------------------------------------------------- Pjump

Source = GranularSourceLinker()
Source:set_particle_density(params.rho_s)
Source:set_fluid_density(params.rho_a)
Source:set_mix_density(Density)
Source:set_gravity(-9.81)
Source:set_cons_gravity(params.consistentRho_in_source)


---------------------------------------------------------------------- Sediment Velocity

Vs = 6.8598478663758
iter = 0
Vs = InterfaceValues:RelVel_ext(params.nu_a*params.rho_a,params.rho_a,params.dp,params.rho_s,9.81,1e-05)
re=InterfaceValues:RE(params.nu_a*params.rho_a,params.rho_a,params.dp,Vs)
Cd =InterfaceValues:CD(re,params.drag_mod)
print(Vs)
print(re)
print(Cd)


DX=1.84*math.pow(1/2,params.numRefs)
RelVel=RelativeVelocityLinker()

RelVel:set_mix_density(Density)
RelVel:set_mix_kinematic_viscosity(Visc)

RelVel:set_vol_limit(params.alpha_min)
RelVel:set_fluid_density(params.rho_a)
RelVel:set_particle_density(params.rho_s)
RelVel:set_particle_diameter(params.dp)
RelVel:set_gravity(-9.81)
RelVel:set_rel_vel(Vs)
RelVel:set_dragCoeff(Cd)
RelVel:activate_relative_vel(false)
RelVel:set_fluid_viscosity(params.nu_a*params.rho_a)
RelVel:set_alpha_max(params.alpha_max)
RelVel:set_phase_parameters(InterfaceValues)
---------------------------------------------------------------------- VelocityGradMag

gamma = ShearStressFV1(approxSpace,u)

---------------------------------------------------------------------- Diffusion

Diff_beta = 0.217
d1= 250e-06
A1= 1.673
k0 = 1.0 + A1 * (1-params.dp/d1)
Diff_factor = Diff_beta*k0*k0*Vs/(2*params.rho_a*9.81)
Diffusion = GranularDiffusionLinker();
Diffusion:set_gamma(gamma)
Diffusion:set_mix_viscosity(Viscosity)
Diffusion:set_diff_factor(Diff_factor)
Diffusion:set_phase_parameters(InterfaceValues)
Diffusion:const_kinematic_visc(true)
print ("	--------------------------Diff	factor	= " .. Diff_factor)


ss_value = math.atan(params.FricMu_2)*180/3.1415926

if params.boolSlipDiff then
	SlipDiff = SlipDiffusion(approxSpace,u)
	SlipDiff:set_theta(ss_value)
	SlipDiff:set_vel(0.05)
	SlipDiff:set_phase_parameters(InterfaceValues)
	if params.boolAveDiff then
		SlipDiff:set_diffusion(Diffusion)
	end
else if params.boolSlipVel then
		SlipVel = SlipVelocity(approxSpace,u)
		SlipVel:set_theta(ss_value)
		SlipVel:set_vel(0.1)
	end
end
------------------------------------------------------------------------------------------
-- Compose the discretization
------------------------------------------------------------------------------------------

-- inner space

NavierStokesDisc = NavierStokesFV1M (fct_cmp_tbl, Inner_total)
NavierStokesDisc:set_exact_jacobian (params.bExactJac)
NavierStokesDisc:set_stokes (params.bStokes)
NavierStokesDisc:set_laplace (params.bNoLaplace)
NavierStokesDisc:set_upwind (params.upwind_m)
NavierStokesDisc:set_upwind_vol(params.upwind_t)
NavierStokesDisc:set_peclet_blend (params.bPecletBlend)
NavierStokesDisc:set_stabilization (params.stab, params.diffLength)


NavierStokesDisc:set_density(Density)
if params.boolRelativeVel then
	NavierStokesDisc:set_relative_velocity(RelVel)
	NavierStokesDisc:set_upwind_rel(params.upwind_r)
end
if params.boolSlipDiff then
	NavierStokesDisc:set_diffusion(SlipDiff)
else
	if params.boolAveDiff then
		NavierStokesDisc:set_diffusion(Diffusion)
	end
	if params.boolSlipVel then
		NavierStokesDisc:set_slip_velocity(SlipVel)
	end
	
end

if (params.boolSource) then
	NavierStokesDisc:set_source(Source)
end
if params.turbViscMethod=="no" then
	NavierStokesDisc:set_kinematic_viscosity (Visc)
else
	NavierStokesDisc:set_kinematic_viscosity(viscosityData)
end
NavierStokesDisc:set_average_gamma(gamma)
NavierStokesDisc:set_phase_parameters(InterfaceValues)

InletDisc = NavierStokesInflowFV1M (NavierStokesDisc)
InletDisc:add ("inflowVel2d","Left,Top")

flowBnd = DirichletBoundary()
flowBnd:add(0.0, "c", "Left")




-- boundary condition at the outlet
OutletDisc = NavierStokesNoNormalStressOutflowFV1M (NavierStokesDisc)
OutletDisc:add ("Right")
OutletDisc:set_phase_parameters(InterfaceValues)

-- boundary condition at the impermeable walls
WallDisc = NavierStokesWall (NavierStokesDisc)
WallDisc:add ("Bottom")




---------------------------------------------------------------------------------------
-- Parameters Inputs
---------------------------------------------------------------------------------------

Density:set_volume_fraction(NavierStokesDisc:volume_fraction())

Visc:set_volume_fraction(NavierStokesDisc:volume_fraction())
Visc:set_eins_viscosity(NavierStokesDisc:einstein_viscosity())
Visc:set_velocity_gradient(NavierStokesDisc:velocity_grad())

if(true) then
    Visc:set_particle_pressure(NavierStokesDisc:particle_pressure())
else
    Visc:set_particle_pressure(NavierStokesDisc:pressure())
end

Diffusion:set_velocity_gradient(NavierStokesDisc:velocity_grad())


RelVel:set_volume_fraction(NavierStokesDisc:volume_fraction())
RelVel:set_volume_grad(NavierStokesDisc:volume_fraction_grad())
RelVel:set_pressure_grad(NavierStokesDisc:pressure_grad())
RelVel:set_einstein_visc(NavierStokesDisc:einstein_viscosity())
RelVel:set_ps_grad(NavierStokesDisc:particle_pressure_grad())



---------------------------------------------------------------------------------------
-- Global Discretization
---------------------------------------------------------------------------------------

-- the global discretization
domainDisc = DomainDiscretization (approxSpace)
domainDisc:add (NavierStokesDisc)
domainDisc:add (InletDisc)
domainDisc:add (OutletDisc)
domainDisc:add (WallDisc)
domainDisc:add(flowBnd)



--domainDisc:add(TransportEq)
--domainDisc:add(OutflowBND)


---------------------------------------------------------------------------------------
-- Time Discretization
---------------------------------------------------------------------------------------

-- create the assembled operator for the solver
local timeDisc = nil
-- create time discretization
if params.timeMethod=="cn" then
	timeDisc = ThetaTimeStep(domainDisc)
	timeDisc:set_theta(0.5) -- Crank-Nicolson method
end
if params.timeMethod=="euler" then
	timeDisc = ThetaTimeStep(domainDisc)
	timeDisc:set_theta(1) -- implicit Euler
end	
if params.timeMethod=="fracstep" then
	timeDisc = ThetaTimeStep(domainDisc,"FracStep")
end
if params.timeMethod=="alex" then
	timeDisc = ThetaTimeStep(domainDisc, "Alexander")
end

------------------------------------------------------------------------------------------
-- Set up the solver
------------------------------------------------------------------------------------------

boolSolution = 1
op, NLSolver, NewtonSolverSteady, limex, boolSolution = myProblem:CreateSolver(domainDisc, approxSpace, timeDisc)



------------------------------------------------------------------------------------------
-- Interpolate initial values
------------------------------------------------------------------------------------------


  --Interpolate(StartValueX, u, "u")
Interpolate(1.0e-5, u, "u")
Interpolate("StartValueY", u, "v")
Interpolate("StartValueP", u, "p")
--Interpolate("StartValueC", u, "c")
Interpolate("ConstValue", u, "c")


------------------------------------------------------------------------------------------
-- Steady State Solution
------------------------------------------------------------------------------------------
time_work_steady=0.0
linsolver_calls = 0
linsolver_steps = 0

if params.doSteadyState and boolSolution == 1 then
	-- Steady state solution.
	time_work_steady, linsolver_calls, linsolver_steps, boolSolution = myProblem:ComputeNonLinearSteadyStateSolution(u, domainDisc, NewtonSolverSteady)
	
	RelVel:activate_relative_vel(params.boolRelativeVel)
    if not(params.boolViscPs) then
        Visc:set_particle_pressure(NavierStokesDisc:pressure())
    end
end

------------------------------------------------------------------------------------------
-- Printing Initial Conditions
------------------------------------------------------------------------------------------
-- start
time = 0
step = 0

	-- write start solution
if boolSolution == 1 then
	print("Writing inittial values")
	out = VTKOutput()
	out:clear_selection()g
	out:select_all(false)
	out:select_nodal ("u,v", "velocity")
	out:select_nodal ("u", "u")
	out:select_nodal ("v", "v")
	out:select_nodal ("p", "p")
	out:select_nodal ("c", "c")
	out:select(Density, "Rho")
	out:select(NavierStokesDisc:einstein_viscosity(), "Mu_eins")
	out:select(Viscosity, "Mu")
	out:select(Viscosity2, "Mu2")
	out:select(RelVel, "RelVel")
	out:select(NavierStokesDisc:particle_pressure(), "Ps")
	out:select(NavierStokesDisc:particle_pressure_grad(), "DPs")
	out:select(gamma, "G")
	out:select(NavierStokesDisc:velocity_grad(), "G2")
	if (params.boolSlipDiff) then
		out:select(SlipDiff, "SDiff")
	else
		if params.boolSlipVel then
			out:select(SlipVel, "SVel")
		end
	end
	out:select(Diffusion, "D")
	out:print_subsets(vtk_file_name, u,allSubsets,step,time)
	print ("Output to file " .. vtk_file_name .. ".vtu  in time t = 0")
	print ("    -   -   -   -   -   -   -   -   -   -   -   -   -   -   ")
	print ("                                                            ")
	print ("    -   -   -   -   -   -   -   -   -   -   -   -   -   -   ")
end


------------------------------------------------------------------------------------------
-- Final Setting
------------------------------------------------------------------------------------------
-- create new grid function for old value
uOld = u:clone()

-- store grid function in vector of  old solutions
solTimeSeries = SolutionTimeSeries()
solTimeSeries:push(uOld, time)


Value_inner1 = Integral(NavierStokesDisc:volume_fraction(), u,"Inner",0.0)
Value_inner2 = Integral(NavierStokesDisc:volume_fraction(), u,"Inner2",0.0)

if (rank == 0 and  boolSolution == 1) then
	myProblem:WriteValues( folder, step, time, Value_inner1, Value_inner2, time_work_steady, 1, 0, linsolver_calls, linsolver_steps,false)
end

NLSolver:add_step_update(viscosityData)
if params.timeMethod == "limex" then
	NLSolver:add_step_update(gamma)

	if params.boolSlipDiff then
		NLSolver:add_step_update(SlipDiff)
	else if params.boolSlipVel then
			NLSolver:add_step_update(SlipVel)
		end
	end
else
	NLSolver:add_inner_step_update(gamma)
	if params.boolSlipDiff then
		NLSolver:add_inner_step_update(SlipDiff)
	else if params.boolSlipVel then
			NLSolver:add_inner_step_update(SlipVel)
		end
	end
	
end



total_Newton_Steps = 0
total_Newton_Steps_fail = 0
total_linsolver_calls_step = 0
total_linsolver_steps_step = 0
tBefore = os.clock()
------------------------------------------------------------------------------------------
-- Time Steps Loop    (Solution)
------------------------------------------------------------------------------------------
if boolSolution == 1 then
	for step = 1, params.numTimeSteps do

		print("++++++ TIMESTEP " .. step .. " BEGIN ++++++")
		tBefore_step = os.clock()
		StartTime = time
		EndTime = time + params.DTmax
		if params.timeMethod == "limex" then
		
			Newton_Steps, Newton_Steps_fail, linsolver_calls_step, linsolver_steps_step, boolSolution  = myProblem:SolveNonlinearProblemLimex(u, limex, NLSolver, StartTime, EndTime)

		else
			Newton_Steps, Newton_Steps_fail, linsolver_calls_step, linsolver_steps_step , dt, boolSolution = myProblem:SolveNonlinearProblem( u, NLSolver, op, timeDisc, solTimeSeries, dt, step,StartTime,EndTime)
		end
		time = EndTime
		tAfter_step = os.clock()
		
		if (step % params.outputFactor == 0 and boolSolution == 1) then
		
			out:print_subsets(vtk_file_name, u,allSubsets,step,time)
			print ("Output to file " .. vtk_file_name .. ".vtu  in time t =  " .. time .. "  Step = " .. step .. ".")
			print(" ")
		end
		if boolSolution == 1 then
			print("++++++ TIMESTEP " .. step .. "  END ++++++")
			print ("    -   -   -   -   -   -   -   -   -   -   -   -   -   -   ")
			print ("                                                            ")
			print ("                                                            ")
			print ("                                                            ")
			print ("    -   -   -   -   -   -   -   -   -   -   -   -   -   -   ")
			print ("                                                            ")
			print ("<<<<<< Total Newton semi   Steps =  " .. Newton_Steps .. "   >>>>>>")
			print ("<<<<<<       Newton success Steps =  " .. Newton_Steps-Newton_Steps_fail .. "   >>>>>>")
			print ("<<<<<<       Newton fail   Steps =  " .. Newton_Steps_fail .. "   >>>>>>")
			print ("                                                            ")
			print ("    -   -   -   -   -   -   -   -   -   -   -   -   -   -   ")
			print ("                                                            ")
			print ("                                                            ")
			print ("                                                            ")
			print ("    -   -   -   -   -   -   -   -   -   -   -   -   -   -   ")
			
			total_Newton_Steps = total_Newton_Steps + Newton_Steps
			total_Newton_Steps_fail = total_Newton_Steps_fail + Newton_Steps_fail
			total_linsolver_calls_step = total_linsolver_calls_step + linsolver_calls_step
			total_linsolver_steps_step = total_linsolver_steps_step + linsolver_steps_step
					
			Value_inner1 = Integral(NavierStokesDisc:volume_fraction(), u,"Inner",0.0)
			Value_inner2 = Integral(NavierStokesDisc:volume_fraction(), u,"Inner2",0.0)
			
			if rank == 0 then
				myProblem:WriteValues( folder, step, time, Value_inner1, Value_inner2, tAfter_step - tBefore_step, Newton_Steps, Newton_Steps_fail, linsolver_calls_step, linsolver_steps_step,false)
			end
		else
			break
		end
		
		
	end
end
tAfter = os.clock()

if boolSolution == 1 then
	------------------------------------------------------------------------------------------
	-- Solution Done
	------------------------------------------------------------------------------------------

	print("-			-")
	print("-------------------------------------------------------------------------------")
	print("Steady state Computation took " .. time_work_steady .. " seconds.")
	print("Temporal Computation took " .. tAfter-tBefore .. " seconds.")
	print("Total Computation took " .. time_work_steady+tAfter-tBefore .. " seconds.")
	print("-------------------------------------------------------------------------------")
	print("")
	print("")
	print ("Output to file " .. vtk_file_name .. ".vtu")
	print("done.")


	if rank == 0 then
		myProblem:WriteValues( folder, params.numTimeSteps, time, Value_inner1, Value_inner2, time_work_steady+tAfter-tBefore, total_Newton_Steps, total_Newton_Steps_fail, total_linsolver_calls_step, total_linsolver_steps_step,true)
	end
end


local Tablename = folder_name .. "/Table_out_" .. numProc ..".csv"
lineWriter = LineWriter()
Headers = " Sim, Vel, H0, W0, Solved\n"
lineWriter:write_line(Tablename,rank_t, Headers, inflow, H_0, W0, boolSolution)

SpaceTimeComm:unsplit()

