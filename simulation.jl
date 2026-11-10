#记得改文件名
using OceanBioME, Oceananigans, Printf, CairoMakie
using Oceananigans.Units

using Random
Random.seed!(11)
grid = RectilinearGrid(size = (128,128, 32),
                       extent = (10kilometer, 10kilometer, 500meters))
coriolis = FPlane(f = 1e-4) # [s⁻¹]
buoyancy = SeawaterBuoyancy()
background_state_parameters = (M = 0,       # s⁻¹, geostrophic shear
                               f = coriolis.f, # s⁻¹, Coriolis parameter
                               N_ = 1e-3,       # s⁻¹, buoyancy frequency
                               H = grid.Lz,
                               g = buoyancy.gravitational_acceleration,
                               α = buoyancy.equation_of_state.thermal_expansion)
T(x, y, z, t, p) = (p.(1/3000)^2 * 5000 + p.N_^2 * (z + p.H)) / (p.g * p.α)
V(x, y, z, t, p) = p.M^2 / p.f * (z + p.H)

V_field = BackgroundField(V, parameters = background_state_parameters)
T_field = BackgroundField(T, parameters = background_state_parameters)
νᵥ = κᵥ = 1e-4 # [m² s⁻¹]
vertical_diffusivity = VerticalScalarDiffusivity(ν = νᵥ, κ = κᵥ)
biogeochemistry = NPZD(grid;
                          scale_negatives = true)

#加surface cooling和bottom drag
Qcool = 0 # W/m² 向上为正
ρ = 1025
cp = 3990
FT = Qcool/(ρ*cp)
T_flux = FluxBoundaryCondition(FT)

bottom_drag = BulkDrag(coefficient =2e-3)
u_bcs = FieldBoundaryConditions(
    bottom = bottom_drag
)
v_bcs = FieldBoundaryConditions(
    bottom = bottom_drag
)

boundary_conditions = (
    T = FieldBoundaryConditions(
        top = T_flux
    ),

    N = FieldBoundaryConditions(),

    P = FieldBoundaryConditions(),

    D = FieldBoundaryConditions(),
    u = u_bcs,
    v = v_bcs

)

model = NonhydrostaticModel(grid;
                            biogeochemistry,
                            advection = WENO(),
                            coriolis,
                            tracers = (:N, :P, :Z,:T,:S ,:D),
                            buoyancy,
                            background_fields = (v = V_field,t=T_field),
                            closure = vertical_diffusivity,

                            boundary_conditions = boundary_conditions
                            )
Ξ(z) = randn() * z / grid.Lz * (z / grid.Lz + 1)

Ũ = 1e-3
uᵢ(x, y, z) = Ũ * Ξ(z)
vᵢ(x, y, z) = Ũ * Ξ(z)

set!(model, u=uᵢ, v=vᵢ, P = 0.03, Z = 0.03, N = 4.0, S = 35, T = 8)
Δx = minimum_xspacing(grid, Center(), Center(), Center())
Δy = minimum_yspacing(grid, Center(), Center(), Center())
Δz = minimum_zspacing(grid, Center(), Center(), Center())
Δt₀ = 0.5 * min(Δx, Δy, Δz) / V(0, 0, 0, 0, background_state_parameters)

simulation = Simulation(model, Δt = Δt₀, stop_time = 15days)
wizard = TimeStepWizard(cfl = 0.5, diffusive_cfl = 0.5, max_Δt = 30minutes)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(5))
progress(sim) = @printf("i: % 6d, sim time: % 10s, wall time: % 10s, Δt: % 10s, CFL: %.2e\n",
                        sim.model.clock.iteration,
                        prettytime(sim.model.clock.time),
                        prettytime(sim.run_wall_time),
                        prettytime(sim.Δt),
                        AdvectiveCFL(sim.Δt)(sim.model))

simulation.callbacks[:progress] = Callback(progress, IterationInterval(20))
u, v, w = model.velocities # unpack velocity `Field`s
ζ = Field(∂x(v) - ∂y(u))
simulation.output_writers[:fields] = JLD2Writer(model, merge(model.tracers, (; u, v, w,ζ));
                                                schedule = TimeInterval(2hours),
                                                filename = "Q_0_m_0_res_128.jld2",
                                                overwrite_existing = true)
run!(simulation)
