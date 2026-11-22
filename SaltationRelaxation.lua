-- SaltationRelaxation.lua
-- Auxiliary scalar convection-reaction problem on the existing distributed
-- Bottom subset, independent of the block-algebra Navier--Stokes problem.
-- Test equation: div(v*q) + q = 1, v = (1,0[,0]).
-- Test inlet: q = 0.2 where x = 0 along Bottom.
-- No changes are made to the main CFD discretization or solver.

local SaltationRelaxation = {}

local function scalar_class(name, suffix)
    local cls = _G[name .. suffix]
    assert(cls, "SaltationRelaxation: missing UG4 class " .. name .. suffix)
    return cls
end

function SaltationRelaxation.Create(dom, dim)
    assert(dom ~= nil, "SaltationRelaxation: domain is required")
    assert(dim == 2 or dim == 3, "SaltationRelaxation: requires 2D or 3D CFD domain")

    print("Creating saltation discretization on Bottom")

    -- Explicit CPU1 DoF layout: the main CFD problem may use CPU6.
    local approxSalt = ApproximationSpace(dom, AlgebraType("CPU", 1))
    local boundarySubsets

	if dim == 2 then
		boundarySubsets = "BottomLeft,BottomRight"
	else
		boundarySubsets = "BottomLeft,BottomRight,BottomFront,BottomBack"
	end

	approxSalt:add_fct( "q", "Lagrange", 1, "Bottom," .. boundarySubsets )
    approxSalt:init_levels()
    approxSalt:init_top_surface()
    approxSalt:print_statistic()

    local suffix = dim .. "dCPU1"
    local q = scalar_class("GridFunction", suffix)(approxSalt)
    q:set(0.0)

    local velocity = dim == 3 and {1.0, 0.0, 0.0} or {1.0, 0.0}
    local saltDisc = ConvectionDiffusionFV1("q", "Bottom")
    saltDisc:set_diffusion(0.0)
    saltDisc:set_velocity(velocity)
    saltDisc:set_reaction_rate(1.0)
    saltDisc:set_source(1.0)
    saltDisc:set_upwind(FullUpwind())

    local domainSalt = scalar_class("DomainDiscretization", suffix)(approxSalt)
    domainSalt:add(saltDisc)
    local saltOutflow = scalar_class("ConvectionDiffusionOutflowFV1", suffix)(saltDisc)
	saltOutflow:add("BottomRight")
	domainSalt:add(saltOutflow)

    -- UG4 conditional Lua Dirichlet callback returns (active, value).
    -- Both endpoints are Bottom vertices in the checked 2D UGX grid.
    local inlet
    if dim == 2 then
        inlet = function(x, y, t)
            return x <= 1.0e-8, 0.2
        end
    else
        inlet = function(x, y, z, t)
            return x <= 1.0e-8, 0.2
        end
    end
    local dirichlet = scalar_class("DirichletBoundary", suffix)()
    dirichlet:add(inlet, "q", "Bottom," .. boundarySubsets)
    domainSalt:add(dirichlet)

    print("Saltation approximation space created")
	return {approxSpace = approxSalt, solution = q, discretization = domainSalt, boundary = dirichlet, boundarySubsets = boundarySubsets, dim = dim}
end

function SaltationRelaxation.Solve(salt)
    local domainSalt = salt.discretization
    local q = salt.solution

    -- AssembledLinearOperator is algebra-specific, not dimension-specific.
    local A = scalar_class("AssembledLinearOperator", "CPU1")(domainSalt)
    local b = q:clone()
    domainSalt:assemble_linear(A, b)
    domainSalt:adjust_solution(q)
    print("SALT_BOTTOM_ASSEMBLY_PASS")

    local solver = scalar_class("BiCGStab", "CPU1")()
    solver:set_preconditioner(scalar_class("ILU", "CPU1")())
    solver:set_convergence_check(scalar_class("ConvCheck", "CPU1")(200, 1.0e-8, 1.0e-8, true))
    assert(solver:init(A, q), "SaltationRelaxation: solver initialization failed")
    assert(solver:apply(q, b), "SaltationRelaxation: linear solve failed")

    print("SALT_BOTTOM_SOLVE_PASS")
    return q
end

function SaltationRelaxation.Validate(salt,vtk_folder)

    local q = salt.solution

    -- Analytical solution
    if salt.dim == 2 then
        _G.SaltationExact = function(x, y, t)
            return 1.0 - 0.8 * math.exp(-x)
        end
    else
        _G.SaltationExact = function(x, y, z, t)
            return 1.0 - 0.8 * math.exp(-x)
        end
    end

    -- Maximum nodal error on Bottom
	local maxError = MaxError("SaltationExact", q, "q", "Bottom," .. salt.boundarySubsets, 0.0)

    print(string.format(
        "SALT_VALIDATION rank=%d max_error=%.12e", ProcRank(), maxError))

    -- Print DoF statistics, including parallel distribution
    salt.approxSpace:print_statistic("component subset interface")

    -- Visualize only the saltation solution
    local vtk = VTKOutput()
    vtk_folder = vtk_folder or "."
	vtk:print_subsets(vtk_folder .. "/SaltationValidation", q, "Bottom," .. salt.boundarySubsets, 0, 0.0)

    print("SALT_VALIDATION_DONE")
end

return SaltationRelaxation
