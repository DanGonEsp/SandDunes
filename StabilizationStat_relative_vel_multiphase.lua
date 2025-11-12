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



geom_name = "quad" -- tri and quad
jumpPressure = false
StatBool = true
boolSource = true
boolGradientPsSource = false
boolAveDiff = true
boolRelativeVel = true


consistentRho_in_source = true



-- Numerical parameters of the discretization
numRefs     = util.GetParamNumber("-numRefs", 1, "number of grid refinements")
numPreRefs     = util.GetParamNumber("-numPreRefs", 0, "number of prerefinements (parallel)")


------------------------------------------------------------------------------------------
-- Get command line parameters
------------------------------------------------------------------------------------------

-- Geometry parameters

if geom_name == "tri" then
	geometry	= util.GetParam ("-geom", "Dune2D_tri_double")
	file_name = "TriDouble"
else
	geometry	= util.GetParam ("-geom", "Dune2D_quads_double")
	file_name = "QuadsDouble"
end
gridName = geometry .. ".ugx"
-- Subsets used in the problem
allSubsets = "Inner, Inner2,Left, Right,Top, Bottom"

------------------------------------------------------------------------------------------
-- Folder and files
------------------------------------------------------------------------------------------

folder = file_name
if not DirectoryExists (folder) then
    CreateDirectory (folder)
end




vtk_file_name = file_name .. "-lev" .. numRefs
if jumpPressure then

    vtk_file_name =file_name .. "-Press"
    value_beta =-0.005         --value_beta = -0.005
    jac = 0.0
    damping_mg = 1.0        --damping_mg = 0.15
else
    vtk_file_name =file_name .. "-NoPress"
    value_beta = -0.01       --0.1--Tri
    jac = 0.0
    damping_mg = 0.5   --0.9  --Tri
end

if boolRelativeVel then
    vtk_file_name =vtk_file_name .. "-RelVel"
else
    vtk_file_name =vtk_file_name .. "-NoRelVel"
end


if  boolGradientPsSource  or boolSource  then

    if boolSource then
        if consistentRho_in_source then vtk_file_name =vtk_file_name .. "-Consistent" end
        vtk_file_name = vtk_file_name .. "-MG_Force"
    end
    
    if boolGradientPsSource then
        vtk_file_name = vtk_file_name .. "-DPs"
    else
        vtk_file_name = vtk_file_name .. "-NoDPs"
    end
    
else
    vtk_file_name = vtk_file_name .. "-NoForce"
end

if boolAveDiff then
    vtk_file_name = vtk_file_name .. "-AveDiff"
else
    vtk_file_name = vtk_file_name .. "-NoAveDiff"
end

if bStokes then
    vtk_file_name = vtk_file_name .. "-Stokes"
end

folder = folder .. "/" .. vtk_file_name
if not DirectoryExists (folder) then
    CreateDirectory (folder)
end

vtk_file_name = folder .. "/" .. vtk_file_name -- VTK output file name base




vtk_file_name = file_name .. "-lev" .. numRefs

if jumpPressure then

    vtk_file_name =file_name .. "-Pres"
    value_beta =-0.005         --value_beta = -0.005
    jac = 0.0
    damping_mg = 1.0        --damping_mg = 0.15
else
    vtk_file_name =file_name .. "-NoPress"
    value_beta = -0.01       --0.1--Tri
    jac = 0.0
    damping_mg = 0.5   --0.9  --Tri
end

if boolRelativeVel then
    vtk_file_name =file_name .. "-RelVel"
end


if  boolGradientPsSource  or boolSource  then

    if boolSource then
        if consistentRho_in_source then vtk_file_name =vtk_file_name .. "-Consistent" end
        vtk_file_name = vtk_file_name .. "-MG_Force"
    end
    
    if boolGradientPsSource then
        vtk_file_name = vtk_file_name .. "-DPs"
    else
        vtk_file_name = vtk_file_name .. "-NoDPs"
    end
    
else
    vtk_file_name = vtk_file_name .. "-NoForce"
end
if boolAveDiff then
    vtk_file_name = vtk_file_name .. "-AveDiff"
end
if bStokes then
    vtk_file_name = vtk_file_name .. "-Stokes"
end



-- Physical parameters
dim 		= util.GetParamNumber("-dim", 2, "dimensionality of the problem")
nu_a 	= util.GetParamNumber("-visc_a", 1.48e-05, "kinematic viscosity")
rho_a 	= util.GetParamNumber("-rho_a", 1.2, "Air Density")
rho_s 	= util.GetParamNumber("-rho_s", 2500, "Sand Density")
dp 	= util.GetParamNumber("-diameter", 1e-03, "Particle Diameter")
nu_s 	= util.GetParamNumber("-visc_s", 1.48e-05, "kinematic viscosity")
inflow		= util.GetParamNumber("-inflow", 10.0, "max. inflow velocity")
c_init		= util.GetParamNumber("-initial concentration", 0.6205, "max volume fraction")


alpha_max		= util.GetParamNumber("-max_concentration", 0.635, "max volume fraction")
alpha_min		= util.GetParamNumber("-min concentration", 0.57, "max volume fraction")
packing_factor        = util.GetParamNumber("-packing_factor", 0.625, "max volume fraction")
viscosity_model		= util.GetParamNumber("-granular_model", 3, "Options:  0 Constant, 1 Proportional, 2 Einstein model, 3 Rheology(I) + Einstein model")
density_model  = util.GetParam("-density_model", "linear", "constant, linear")
interface_value  = util.GetParamNumber("-interface_value",  10, "interface value")
FR = 0.05
B_phi = 1
deltaGamma = 1e-07;
Visc_limit = 1e15

deltaPs = 1.48e-04;
deltaI = 1e-03;
FricMu_1=0.38
FricMu_2=0.64
I_0 = 0.279

drag_mod = util.GetParamNumber("-drag_model", 2, "Options:  0 StokesLaw, 1 formula, 2 Schiller-Naumann, 3 Turton and Levenspiel") --Model 0 pow(0.63+4.8/sqrt(RE),2.0);

    


bStokes 	= util.GetParamNumber("-Stokes", false ,"If defined, only Stokes Eq. computed")
bNoLaplace 	= util.GetParamNumber("-noLaplace", false,"If defined, only laplace term used")
bExactJac 	= util.GetParamNumber("-exactJac", jac,"If defined, exact jacobian used")
bPecletBlend= util.GetParamNumber("-PecletBlend", false,"If defined, Peclet Blend used")
upwind_m      = util.GetParam("-upwind_m", "full", "Upwind type full or lps")
upwind_t      = util.GetParam("-upwind_t", "full", "Upwind type full or lps")
bPac        = util.GetParamNumber("-pac", false,"If defined, pac upwind used")
diffLength  = util.GetParam("-difflength", "raw", "fivepoint, raw, corDiffusion length type")
stab        = util.GetParam("-stab", "fields_2", "Stabilization type (fields or flow viscosity or karimian)")

time_years = util.GetParamNumber("-dT_y",0.0)
time_days = util.GetParamNumber("-dT_d",0.0)
time_hours = util.GetParamNumber("-dT_h",0.0)
time_seconds = util.GetParamNumber("-dT_ss",2.0)

dt_s =  31536000 * time_years + 86400*time_days +  3600*time_hours+ time_seconds
DTmax= util.GetParamNumber("-DTmax", dt_s, "max  DT")
DTmin= util.GetParamNumber("-DTmin", 0.001, "min  DT")
dt = util.GetParamNumber("-dt",  1.0)

CFL_max= util.GetParamNumber("-cfl", 100, "max  CFL number")
UpdateDt     = util.GetParamNumber("-UpdateDT", 5)
modifyDT     = util.GetParamNumber("-modifyDT", true)

incr_factor     = util.GetParamNumber("-CFL_factor", 1.15)
red_factor_fail     = util.GetParamNumber("-CFL_factor", 0.7)
red_factor_success     = util.GetParamNumber("-CFL_factor", 0.8)

maxConvRate = util.GetParamNumber("-maxConvRate", 0.85)
minConvRate = util.GetParamNumber("-minConvRate", 0.7)

numTimeSteps =  util.GetParamNumber("-numTimeSteps", 20	)
outputFactor     = util.GetParam("-output", 1, "output every ... steps")




-- Parameters of the solver

max_newton_steps_steady_state=util.GetParamNumber("-numNewtonSteps", 100)
max_newton_steps_transcient=util.GetParamNumber("-numNewtonSteps", 50)
max_linear_steps=util.GetParamNumber("-numLinearIter", 200)
timeMethod = util.GetParam("-timeMethod","euler","cn euler   fracstep   alex")




turbViscMethod = util.GetParam("-turbulenceModel","no","TurbVismodel type no , dyn or sma")
modellconstant = util.GetParamNumber("-c",0.1)







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
print ("	upwind_m	= " .. upwind_m)
print ("	upwind_t	= " .. upwind_t)
print ("	pac			= " .. tostring (bPac))
print ("	stab		= " .. stab)
print ("	difflength	= " .. diffLength)

------------------------------------------------------------------------------------------
-- Initialize UG4, load, refine and distribute the grid
------------------------------------------------------------------------------------------

InitUG (dim, AlgebraType("CPU", dim +2))

if dim == 3 then
	fct_cmp_tbl = {"u", "v", "w", "p", "c"}
	vel_cmp_tbl = {"u", "v", "w"}
else
	fct_cmp_tbl = {"u", "v", "p", "c"}
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



--[[function VolumeFraction2(x,y)
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

function VolumeFraction2(x,y)
	return  c_init
end

function InitialValue_FractionVolume(x,y)
	--if   y < 5 and x>-1 and x<50  then 
		value= VolumeFraction(x,y)
	--else
		--value = VolumeFraction2(x,y)
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

---------------------------------------------------------------------- Interface Properties

InterfaceValues = Interface()

InterfaceValues:set_particle_diameter(dp)
InterfaceValues:set_particle_density(rho_s)
InterfaceValues:set_air_density(rho_a)
InterfaceValues:set_fluid_Visc(nu_a*rho_a)
InterfaceValues:set_alpha_max(alpha_max)
InterfaceValues:set_alpha_min(alpha_min)
InterfaceValues:set_FR(FR)
InterfaceValues:set_B_phi(B_phi)
InterfaceValues:set_deltaGamma(deltaGamma)
InterfaceValues:set_limit(Visc_limit)
InterfaceValues:set_bool_particle_pressure_force(boolGradientPsSource)
InterfaceValues:set_bool_consistent_gravity(false)
InterfaceValues:set_time_step_factor(dt)
InterfaceValues:set_drag_model(drag_mod)
InterfaceValues:set_bool_initialized(true)

---------------------------------------------------------------------- Density

Density = GranularDensityLinker(); 
Density:set_fluid_density(rho_a)
Density:set_particle_density(rho_s)
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
--Visc:set_packing_factor(packing_factor)
Visc:set_mix_density(Density)
Visc:set_interface_volume_fraction(interface_value)
Visc:set_limit(Visc_limit)
Visc:set_deltaPs(deltaPs)
Visc:set_deltaI(deltaI)
Visc:set_FricMu_1(FricMu_1)
Visc:set_FricMu_2(FricMu_2)
Visc:set_I_0(I_0)
Visc:set_deltaGamma(deltaGamma)



if turbViscMethod=="dyn" then
	viscosityData = FV1DynamicTurbViscData(approxSpace,u)
	
else 
	viscosityData = FV1SmagorinskyTurbViscData(approxSpace,u,modellconstant)
 end
viscosityData:set_turbulence_zero_bnd("Left")
viscosityData:set_kinematic_viscosity(Visc)


---------------------------------------------------------------------- Pjump
PjumpShape= JumpShapeLinker()
PjumpShape:set_interface_volume_fraction(interface_value)


normal = InterfaceNormalLinker()
normal:set_interface_volume_fraction(interface_value)

--GradientPsSource = ParticlePressureGradientLinker()
--GradientPsSource:set_interface_volume_fraction(interface_value)


Source = GranularSourceLinker()
Source:set_particle_density(rho_s)
Source:set_fluid_density(rho_a)
Source:set_mix_density(Density)
Source:set_gravity(-9.81)
--Source:set_packing_factor(packing_factor)
Source:set_cons_gravity(consistentRho_in_source)


---------------------------------------------------------------------- Sediment Velocity

Ws = 6.8598478663758
iter = 0
Ws = InterfaceValues:RelVel_ext(nu_a*rho_a,rho_a,dp,rho_s,9.81,1e-05)
re=InterfaceValues:RE(nu_a*rho_a,rho_a,dp,Ws)
Cd =InterfaceValues:CD(re,drag_mod)
print(Ws)
print(re)
print(Cd)


DX=1.84*math.pow(1/2,numRefs)
W=RelativeVelocityLinker()

W:set_mix_density(Density)
W:set_mix_kinematic_viscosity(Visc)
W:set_gravity_force(Source)

W:set_vol_limit(alpha_min)
W:set_fluid_density(rho_a)
W:set_particle_density(rho_s)
W:set_particle_diameter(dp)
W:set_gravity(-9.81)
W:set_rel_vel(Ws)
W:set_dragCoeff(Cd)
W:activate_relative_vel(boolRelativeVel)
W:set_fluid_viscosity(nu_a*rho_a)
W:set_alpha_max(alpha_max)
W:set_phase_parameters(InterfaceValues)
---------------------------------------------------------------------- VelocityGradMag



gamma = ShearStressFV1(approxSpace,u)



---------------------------------------------------------------------- Diffusion

Diff_beta = 0.217
d1= 250e-06
A1= 1.673
k0 = 1.0 + A1 * (1-dp/d1)
Diff_factor = Diff_beta*k0*k0*Ws*nu_a/9.81
Diffusion = GranularDiffusionLinker(); 
Diffusion:set_gamma(gamma)
Diffusion:set_diff_factor(Diff_factor)
print ("	--------------------------Diff	factor	= " .. Diff_factor)

------------------------------------------------------------------------------------------
-- Compose the discretization
------------------------------------------------------------------------------------------

-- inner space

NavierStokesDisc = NavierStokesFV1M (fct_cmp_tbl, {"Inner","Inner2"})
NavierStokesDisc:set_exact_jacobian (bExactJac)
NavierStokesDisc:set_stokes (bStokes)
NavierStokesDisc:set_laplace ( bNoLaplace)
NavierStokesDisc:set_upwind (upwind_m)
NavierStokesDisc:set_upwind_vol(upwind_t) 
NavierStokesDisc:set_peclet_blend (bPecletBlend)
NavierStokesDisc:set_stabilization (stab, diffLength)
--NavierStokesDisc:set_pac_upwind (bPac)


NavierStokesDisc:set_density(Density)
NavierStokesDisc:set_relative_velocity(W)
NavierStokesDisc:set_interface_value(interface_value)
NavierStokesDisc:set_phase_parameters(InterfaceValues)
if boolAveDiff then
	NavierStokesDisc:set_diffusion(Diffusion)
end
if (boolSource) then
	NavierStokesDisc:set_source(Source)
end

if turbViscMethod=="no" then
	NavierStokesDisc:set_kinematic_viscosity (Visc)
else
	NavierStokesDisc:set_kinematic_viscosity(viscosityData)
end



InletDisc = NavierStokesInflowFV1M (NavierStokesDisc)
--InletDisc:add ("inflowVel3d", "Inlet,Top,Bottom")
InletDisc:add ("inflowVel2d","Left,Top")

-- boundary condition at the outlet
OutletDisc = NavierStokesNoNormalStressOutflowFV1M (NavierStokesDisc)
OutletDisc:add ("Right")
OutletDisc:set_phase_parameters(InterfaceValues)

-- boundary condition at the impermeable walls
WallDisc = NavierStokesWall (NavierStokesDisc)
WallDisc:add ("Bottom")

--PressureOutlet = DirichletBoundary()
--PressureOutlet:add("PressureBoundary", "p", "Bottom")



--Stress = NavierStokesInflowStressFV1(NavierStokesDisc)
--Stress:add("Left,Top")
--Stress:set_velocity(Vel)



---------------------------------------------------------------------------------------
-- Viscosity, Density and Gravitation Input
---------------------------------------------------------------------------------------

Density:set_volume_fraction(NavierStokesDisc:volume_fraction())

Visc:set_volume_fraction(NavierStokesDisc:volume_fraction())
Visc:set_mix_viscosity(NavierStokesDisc:einstein_viscosity())
Visc:set_velocity_gradient(NavierStokesDisc:velocity_grad())
Visc:set_particle_pressure(NavierStokesDisc:particle_pressure())


PjumpShape:set_volume_fraction(NavierStokesDisc:volume_fraction())



if (jumpPressure) then
	NavierStokesDisc:set_jump_shape(PjumpShape,diffLength)
	NavierStokesDisc:set_interface_normal(normal)
	--NavierStokesDisc:set_vol_fraction(VolFraction)
end

normal:set_volume_fraction(NavierStokesDisc:volume_fraction())
normal:set_volume_grad(NavierStokesDisc:volume_fraction_grad())


W:set_volume_fraction(NavierStokesDisc:volume_fraction())
W:set_volume_grad(NavierStokesDisc:volume_fraction_grad())
W:set_pressure_grad(NavierStokesDisc:pressure_grad())
W:set_einstein_visc(NavierStokesDisc:einstein_viscosity())
W:set_ps_grad(NavierStokesDisc:particle_pressure_grad())

--GradientPsSource:set_particle_pressure(NavierStokesDisc:particle_pressure())

---------------------------------------------------------------------- Scale

Viscosity = ScaleLinker();
Viscosity:set_import_1(Density)
Viscosity:set_import_2(Visc)

Viscosity2 = ScaleLinker();
Viscosity2:set_import_1(Density)
Viscosity2:set_import_2(viscosityData)



---------------------------------------------------------------------- Scale

-- the global discretization
domainDisc = DomainDiscretization (approxSpace)
domainDisc:add (NavierStokesDisc)
domainDisc:add (InletDisc)
domainDisc:add (OutletDisc)
domainDisc:add (WallDisc)
--domainDisc:add (PressureOutlet)



--domainDisc:add(TransportEq)
--domainDisc:add(OutflowBND)




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
			rim = true,
			cycle = "V",
			smoother =
			{
				type = "ilu",
				beta = value_beta,
				damping 	= damping_mg,
				sort	= false,
				--sortEps 	= 1.e-50,
				inversionEps 	= 1.e-16,
				--consistentInterfaces   = false
				--overlap 		= false,
				--ordering 		= nil
			},
			preSmooth = 3,
			postSmooth = 3,
			baseSolver = "lu",
			baseLevel = numPreRefs
		},
		convCheck =
		{
			type		= "standard",
			iterations	= max_linear_steps,
			absolute	= 1e-12,
			reduction	= 1e-3,
			verbose		= true
		}
	},
	lineSearch =
	{
		type			= "standard",
		maxSteps		=12,
		lambdaStart		= 1.5625,
		lambdaReduce	= 0.8,
		acceptBest 		= true,
		checkAll		= false,
		--suffDesc		= 0.25,
		maxDefect	= 2e20
		
	},
	convCheck =
	{
		type		= "standard",
		iterations	= max_newton_steps_steady_state,
		absolute	= 1e-6,
		reduction	= 1e-7,
		verbose		= true
	}
}

solver = util.solver.CreateSolver(solverDesc)

	
------------------------------------------------------------------------------------------
-- Prepare the initial guess for the pressure
------------------------------------------------------------------------------------------
if StatBool then
	--Interpolate(StartValueX, u, "u")
	Interpolate(1.0e-5, u, "u")
	Interpolate("StartValueY", u, "v")
	Interpolate("StartValueP", u, "p")
	Interpolate("InitialValue_FractionVolume", u, "c")

	-- grid function for the solution

	-- Fix the mass fraction and solve the linear problem for the momentum

	fixer = DirichletBoundary()
	domainDisc:add(fixer)
	fixer:invert_subset_selection()
	fixer:add("c", "")

	solver:init(AssembledOperator(domainDisc))
	solver:add_inner_step_update(viscosityData)

	
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
    num_linsolver_calls = solver:num_linsolver_calls(1)
    num_linsolver_steps = solver:num_linsolver_steps(1)
    average_linear_steps = solver:average_linear_steps(1)
    total_linsolver_calls = solver:total_linsolver_calls()
    total_linsolver_steps = solver:total_linsolver_steps()
    total_average_linear_steps = solver:total_average_linear_steps()
    last_num_newton_steps = solver:last_num_newton_steps()
    total_average_non_linear_rates = solver:total_average_non_linear_rates()
    
    print("num_newton_steps = " .. num_newton_steps .. ".")
    print("num_linsolver_calls = " .. num_linsolver_calls .. ".")
    print("num_linsolver_steps = " .. num_linsolver_steps .. ".")
    print("average_linear_steps = " .. average_linear_steps .. ".")
    print("total_linsolver_calls = " .. total_linsolver_calls .. ".")
    print("total_linsolver_steps = " .. total_linsolver_steps .. ".")
    print("total_average_linear_steps = " .. total_average_linear_steps .. ".")
    print("last_num_newton_steps = " .. last_num_newton_steps .. ".")
    print("total_average_non_linear_rates = " .. total_average_non_linear_rates .. ".")
    
    --solver:print_average_convergence()
    solver:clear_average_convergence();
    
    
	print("Computation for steady state took " .. tAfter_s-tBefore_s .. " seconds.")
	domainDisc:remove (fixer)
	
	print("++++++++++++++++++++++++ INITIAL CONDITIONS  (STEADY STATE DONE) ++++++++++++++++++++++++")	
	
	
	
	
	W:activate_relative_vel(boolRelativeVel)
	
else

	Interpolate(0.0, u, "u")
	Interpolate("StartValueY", u, "v")
	Interpolate("StartValueP", u, "p")
	Interpolate("InitialValue_FractionVolume", u, "c")
	
end

convCheck =
{
	type		= "standard",
	iterations	= max_newton_steps_transcient,
	absolute	= DTmax*1e-6,
	reduction	= 1e-7,
	verbose		= true
}

solverDesc.convCheck = convCheck
solver = util.solver.CreateSolver(solverDesc)


------------------------------------------------------------------------------------------
-- Apply the solver
------------------------------------------------------------------------------------------
-- start
time = 0
step = 0

	-- write start solution
print("Writing inittial values")
out = VTKOutput()
out:clear_selection()
out:select_all(false)
out:select_nodal ("u,v", "velocity")
out:select_nodal ("u", "u")
out:select_nodal ("v", "v")
out:select_nodal ("p", "p")
out:select_nodal ("c", "c")
out:select(Density, "Rho")
out:select(Viscosity, "Mu")
out:select(Viscosity2, "Mu2")
out:select(W, "W")
out:select(NavierStokesDisc:particle_pressure(), "Ps")
out:select(NavierStokesDisc:particle_pressure_grad(), "DPs")
out:select(Diffusion, "D")
out:select(gamma, "G")
out:print_subsets(vtk_file_name, u,allSubsets,step,time)
print ("Output to file " .. vtk_file_name .. ".vtu  in time t = 0")
print ("    -   -   -   -   -   -   -   -   -   -   -   -   -   -   ")
print ("                                                            ")
print ("    -   -   -   -   -   -   -   -   -   -   -   -   -   -   ")
print ("                                                            ")
print ("    -   -   -   -   -   -   -   -   -   -   -   -   -   -   ")
print ("                                                            ")
print ("    -   -   -   -   -   -   -   -   -   -   -   -   -   -   ")
print ("                                                            ")
print ("    -   -   -   -   -   -   -   -   -   -   -   -   -   -   ")


solver:init(op)
solver:add_step_update(viscosityData)

if boolAveDiff then
	solver:add_inner_step_update(gamma)
end


if solver:prepare(u) == false then
	print ("Newton solver prepare failed.") exit()
end
    
-- create new grid function for old value
uOld = u:clone()

tBefore = os.clock()

-- store grid function in vector of  old solutions
solTimeSeries = SolutionTimeSeries()
solTimeSeries:push(uOld, time)


local file = io.open(folder .. "/Newton_Iterations.txt", "w+")
file:write("Step" .. " \t " .. "Time" .. " \t " .. "TNSteps" .. " \t " .. "SNSteps" .. " \t " .. "FNSteps" .. " \n")
file:write(" " .. step .. " \t " .. time .. " \t " .. 1 .. " \t " .. 1 .. " \t " .. 0 .. " \n")
file:close()

Value_inner1 = Integral(NavierStokesDisc:volume_fraction(), u,"Inner",0.0)
Value_inner2 = Integral(NavierStokesDisc:volume_fraction(), u,"Inner2",0.0)


local file = io.open(folder .. "/Integral.txt", "w+")
file:write("Step" .. " \t " .. "Time" .. " \t " .. "Vol-Dom_1" .. " \t " .. "Vol-Dom_1" .. " \n")
file:write(" " .. step .. " \t " .. time .. " \t " .. Value_inner1 .. " \t " .. Value_inner2 .. " \n")
file:close()

for step = 1, numTimeSteps do
	print("++++++ TIMESTEP " .. step .. " BEGIN ++++++")
	
    Newton_Steps = 0
    Newton_Steps_fail = 0
    CompletedStep = false
    time2 = time
	while CompletedStep==false  do
        Newton_Steps = Newton_Steps+1
		-- choose time step

		do_dt = math.min(dt,math.max((time+DTmax-time2), 0.0))
  
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
            dt = math.max(dt*red_factor_fail,0.99999*DTmin)
            print("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<")
            print("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<   Reducing Timestep in step  " .. step-1 + (time2+do_dt-time)/DTmax .. ".")
            print("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<   new DT          =    " .. dt .. "     Time = " .. time2 .."")
            print("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<   DTmax       =    " .. DTmax .. "     DTmin = " .. DTmin .."")
            print("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<   Rho         =    " .. solver:total_average_non_linear_rates() .. ".")
            print("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<   Dt_factor   =    " .. red_factor_fail .. ".")
            
			
			if dt < DTmin  or modifyDT== false then
				print ("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx            Time step below minimum. Aborting. Failed at point with step  " .. step-1+ (time2+do_dt-time)/DTmax .. "."); exit();

			else
				VecScaleAssign(u, 1.0, solTimeSeries:latest())
			
			end
            --solver:print_average_convergence()
            solver:clear_average_convergence();
		else
            
            if  (time2 + do_dt + (1e-07) * DTmin -time)/DTmax > 1.0 then
            
                time= timeDisc:future_time()                                        -- update new time
                                                                                
                oldestSol = solTimeSeries:oldest()                                  -- get oldest solution
                                                                                    
                VecAssign(oldestSol, u)                                             -- copy values into oldest solution (we reuse the memory here)
                                                                                    
                solTimeSeries:push_discard_oldest(oldestSol, time)                  -- push oldest solutions with new values to front, oldest sol pointer is poped from end
                
                CompletedStep = true
                
                
            else
                
                time2 = timeDisc:future_time()                                      -- update new time
                                                                                    
                oldestSol = solTimeSeries:oldest()                                  -- get oldest solution
                                                                                    
                VecAssign(oldestSol, u)                                             -- copy values into oldest solution (we reuse the memory here)
                                                                                    
                solTimeSeries:push_discard_oldest(oldestSol, time2)                 -- push oldest solutions with new values to front, oldest sol pointer is poped from end
                
                CompletedStep = false
                
            
            end
            
            total_average_non_linear_rates = solver:total_average_non_linear_rates()

            --solver:print_average_convergence()
            solver:clear_average_convergence();
            
            if CompletedStep then
                frac_step = 0
                Local_Time = time
                Dt_factor = dt/do_dt
            else
                frac_step = -1 + (time2-time)/DTmax
                Local_Time = time2
            end
            
            
            if modifyDT then
                if(total_average_non_linear_rates>maxConvRate) then
                    dt = math.max(do_dt*red_factor_success,1.00001*DTmin)
                    print ("-------------------------------------------------------------------Time step decrease at Step " .. step-1 + (time2-time)/DTmax .. ", dt =  " .. dt .. ". ")
                else if CompletedStep== false then
                    if(total_average_non_linear_rates<minConvRate ) then
                        dt=math.min(incr_factor*dt,DTmax)
                        print ("-------------------------------------------------------------------Time step increased at Step " .. step-1 + (time2-time)/DTmax ..", dt =  " .. dt .. ". ")
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
            print("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<   new DT    =    " .. dt .. "     Time = " .. Local_Time .."")
            print("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<   old DT    =    " .. do_dt .."")
            print("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<   DTmax =    " .. DTmax .. "     DTmin = " .. DTmin .."")
            print("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<   Dt_factor   =    " .. dt/do_dt .. ".")
            print("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<   Rho   =    " .. total_average_non_linear_rates .. ".")
                        

		end
			
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
		out:select(Viscosity, "Mu")
		out:select(gamma, "G")
		out:select(Viscosity2, "Mu2")
		out:select(W, "W")
		out:select(NavierStokesDisc:particle_pressure(), "Ps")
		out:select(NavierStokesDisc:particle_pressure_grad(), "DPs")
		out:select(Diffusion, "D")

		

		out:print_subsets(vtk_file_name, u,allSubsets,step,time)
		print ("Output to file " .. vtk_file_name .. ".vtu  in time t =  " .. time .. "  Step = " .. step .. ".")
		print(" ")
	end
 
 
        -- Save number of newton iterations for every time step
    file = io.open(folder .. "/Newton_Iterations.txt", "a")
    file:write(" " .. step .. " \t " .. time .. " \t " .. Newton_Steps .. " \t " .. Newton_Steps-Newton_Steps_fail .. " \t " .. Newton_Steps_fail .. " \n")
    file:close()
    
    Value_inner1 = Integral(NavierStokesDisc:volume_fraction(), u,"Inner",0.0)
    Value_inner2 = Integral(NavierStokesDisc:volume_fraction(), u,"Inner2",0.0)

    file = io.open(folder .. "/Integral.txt", "a")
    file:write(" " .. step .. " \t " .. time .. " \t " .. Value_inner1 .. " \t " .. Value_inner2 .. " \n")
    file:close()

    
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



