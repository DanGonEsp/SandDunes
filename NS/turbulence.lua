
--------------------------------------------------------------------------------
-- RANS turbulence FV1 test
--
-- Current model:
--     dk/dt     = 0
--     domega/dt = 0
--
-- Expected result:
--     k(x,t)     = k(x,0)
--     omega(x,t) = omega(x,0)
--------------------------------------------------------------------------------

ug_load_script("ug_util.lua")
ug_load_script("util/domain_disc_util.lua")

RequiredPlugins({"NavierStokes"})

--------------------------------------------------------------------------------
-- Parameters
--------------------------------------------------------------------------------

dim = util.GetParamNumber("-dim", 2, "problem dimension")

gridName = util.GetParam("-geom", "cylinderp.ugx")

numRefs = util.GetParamNumber("-numRefs", 2, "number of refinements")

startTime = util.GetParamNumber("-start", 0.0)
endTime = util.GetParamNumber("-end", 10.0)

numTimeSteps = util.GetParamNumber("-numTimeSteps", 10)

dt = (endTime - startTime) / numTimeSteps

outputFactor = util.GetParamNumber("-output", 1)
--------------------------------------------------------------------------------
--Folder
--------------------------------------------------------------------------------
folder = "RANSTurbulenceTest"

if ProcRank() == 0 then
    if not DirectoryExists(folder) then
        CreateDirectory(folder)
    end
end

SynchronizeProcesses()

vtk_file_name = folder .. "/RANSTurbulence"

--------------------------------------------------------------------------------
-- Subsets
--------------------------------------------------------------------------------

allSubsets = "Inner,Inlet,Outlet,UpperWall,LowerWall,CylinderWall"

--------------------------------------------------------------------------------
-- Initialize UG4
--------------------------------------------------------------------------------

-- Only two unknowns:
--     k
--     omega
InitUG(dim, AlgebraType("CPU", 1))

--------------------------------------------------------------------------------
-- Domain
--------------------------------------------------------------------------------

dom = util.CreateDomain(gridName, 0)

for i = 1, numRefs do
    util.refinement.CreateRegularHierarchy(dom, 1, true)
end

print("Domain info:")
print(dom:domain_info():to_string())

--------------------------------------------------------------------------------
-- Approximation space
--------------------------------------------------------------------------------

approxSpace = ApproximationSpace(dom)

approxSpace:add_fct("k", "Lagrange", 1, allSubsets)
approxSpace:add_fct("omega", "Lagrange", 1, allSubsets)

approxSpace:init_levels()
approxSpace:init_top_surface()
approxSpace:print_statistic()

util.solver.defaults.approxSpace = approxSpace

--------------------------------------------------------------------------------
-- Grid function
--------------------------------------------------------------------------------

u = GridFunction(approxSpace)
u:set(0.0)

--------------------------------------------------------------------------------
-- RANS turbulence discretization
--------------------------------------------------------------------------------
velGrad = ConstUserMatrix2d(0.0)
velGrad:set_entry(0, 0, 10.0)
wallDistance = ConstUserNumber(1.01)
nu = ConstUserNumber(1.48e-05)
TurbulenceDisc = RANSTurbulenceFV1("k,omega", "Inner")
TurbulenceDisc:set_velocity({0.0, 0.0})
TurbulenceDisc:set_upwind("full")
TurbulenceDisc:set_velocity_gradient(velGrad)
TurbulenceDisc:set_wall_distance(wallDistance)
TurbulenceDisc:set_kinematic_viscosity(nu)
nuT = TurbulenceDisc:turbulent_kinematic_viscosity()
nuEff = ScaleAddLinkerNumber()
nuEff:add(nuT,1.0)
nuEff:add(nu,1.0)


domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(TurbulenceDisc)

--------------------------------------------------------------------------------
-- Boundary conditions
--
-- For this first mass-only test, keep k and omega fixed on the outer boundary.
--------------------------------------------------------------------------------

bnd = DirichletBoundary()

bnd:add(2.0, "k", "Inlet,Outlet,UpperWall,LowerWall,CylinderWall")
bnd:add(3.0, "omega", "Inlet,Outlet,UpperWall,LowerWall,CylinderWall")

domainDisc:add(bnd)

--------------------------------------------------------------------------------
-- Time discretization
--------------------------------------------------------------------------------

timeDisc = ThetaTimeStep(domainDisc)
timeDisc:set_theta(1.0)

--------------------------------------------------------------------------------
-- Solver
--------------------------------------------------------------------------------

linSolver = BiCGStab()

precond = ILU()
linSolver:set_preconditioner(precond)

convCheck = ConvCheck()
convCheck:set_maximum_steps(200)
convCheck:set_minimum_defect(1e-12)
convCheck:set_reduction(1e-10)
convCheck:set_verbose(true)

linSolver:set_convergence_check(convCheck)

newton = NewtonSolver()
newton:set_linear_solver(linSolver)

lineSearch = StandardLineSearch(5, 1.0, 0.5, true)
newton:set_line_search(lineSearch)


newtonConvCheck = ConvCheck(100, 1e-10, 1e-8, true)
newton:set_convergence_check(newtonConvCheck)

op = AssembledOperator(timeDisc)

newton:init(op)

--------------------------------------------------------------------------------
-- Initial conditions
--------------------------------------------------------------------------------
function StartK(x, y)
    return 1.0 + x
end

function StartOmega(x, y)
    return 2.0 + x
end

Interpolate("StartK", u, "k")
Interpolate("StartOmega", u, "omega")

--------------------------------------------------------------------------------
-- Output
--------------------------------------------------------------------------------

out = VTKOutput()

out:clear_selection()
out:select_all(false)

out:select_nodal("k", "k")
out:select_nodal("omega", "omega")
out:select(nuT, "nu_t")
out:select(nuEff, "nuEff")



time = startTime
step = 0

out:print_subsets(vtk_file_name, u,allSubsets,step,time, true)

--------------------------------------------------------------------------------
-- Solution time series
--------------------------------------------------------------------------------

uOld = u:clone()

solTimeSeries = SolutionTimeSeries()
solTimeSeries:push(uOld, time)

--------------------------------------------------------------------------------
-- Time loop
--------------------------------------------------------------------------------

for step = 1, numTimeSteps do

    print("")
    print("------------------------------------------------------------")
    print("TIMESTEP " .. step)
    print("------------------------------------------------------------")

    timeDisc:prepare_step(solTimeSeries, dt)

    newton:prepare(u)

    if newton:apply(u) == false then
        print("Newton solver failed.")
        exit()
    end

    time = time + dt

    uOld = u:clone()
    solTimeSeries:push_discard_oldest(uOld, time)

    if step % outputFactor == 0 then
       out:print_subsets(vtk_file_name, u,allSubsets,step,time, true)
    end
end

out:write_time_pvd(vtk_file_name, u)

print("")
print("RANS turbulence mass-term test completed.")