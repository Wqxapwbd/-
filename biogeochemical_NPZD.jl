using OceanBioME, Oceananigans, Printf, CairoMakie
using Oceananigans.Units
using Random
Random.seed!(11)

grid = RectilinearGrid(size = (32, 32, 8),
                       extent = (100kilometer, 100kilometer, 500meters))
coriolis = FPlane(f = 1e-4) # [s⁻¹]
buoyancy = SeawaterBuoyancy()
background_state_parameters = (M = 1e-4,       # s⁻¹, geostrophic shear
                               f = coriolis.f, # s⁻¹, Coriolis parameter
                               N_ = 1e-4,       # s⁻¹, buoyancy frequency
                               H = grid.Lz,
                               g = buoyancy.gravitational_acceleration,
                               α = buoyancy.equation_of_state.thermal_expansion)
T(x, y, z, t, p) = (p.M^2 * x + p.N_^2 * (z + p.H)) / (p.g * p.α)
V(x, y, z, t, p) = p.M^2 / p.f * (z + p.H)

V_field = BackgroundField(V, parameters = background_state_parameters)
T_field = BackgroundField(T, parameters = background_state_parameters)
νᵥ = κᵥ = 1e-4 # [m² s⁻¹]
vertical_diffusivity = VerticalScalarDiffusivity(ν = νᵥ, κ = κᵥ)
biogeochemistry = NPZD(grid)


model = NonhydrostaticModel(grid;
                            biogeochemistry,
                            advection = WENO(),
                            tracers = (:N, :P, :Z,:T,:S ,:D),
                            timestepper = :RungeKutta3,
                            coriolis,
                            buoyancy,
                            background_fields = (T = T_field, v = V_field),
                            closure = vertical_diffusivity)
Ξ(z) = randn() * z / grid.Lz * (z / grid.Lz + 1)
N=model.tracers.N
Z = model.tracers.Z
P=model.tracers.P
D=model.tracers.D




Ũ = 1e-5
uᵢ(x, y, z) = Ũ * Ξ(z)+ 1e-6 * randn()
vᵢ(x, y, z) = Ũ * Ξ(z)+ 1e-6 * randn()

set!(model, u=uᵢ, v=vᵢ, P = 0.5, Z = 0.05, N= 3.5,D=0.1, S = 35, T =20)
Δx = minimum_xspacing(grid, Center(), Center(), Center())
Δy = minimum_yspacing(grid, Center(), Center(), Center())
Δz = minimum_zspacing(grid, Center(), Center(), Center())

Δt₀ = 0.3 * min(Δx, Δy, Δz) / V(0, 0, 0, 0, background_state_parameters)

simulation = Simulation(model, Δt = Δt₀, stop_time = 3days)
wizard = TimeStepWizard(cfl = 0.3, diffusive_cfl = 0.3, max_Δt = 20minutes)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))
progress(sim) = @printf("i: % 6d, sim time: % 10s, wall time: % 10s, Δt: % 10s, CFL: %.2e\n",
                        sim.model.clock.iteration,
                        prettytime(sim.model.clock.time),
                        prettytime(sim.run_wall_time),
                        prettytime(sim.Δt),
                        AdvectiveCFL(sim.Δt)(sim.model))

simulation.callbacks[:progress] = Callback(progress, IterationInterval(20))
u, v, w = model.velocities # unpack velocity `Field`s
ζ = Field(∂x(v) - ∂y(u))

simulation.output_writers[:fields] = JLD2Writer(model,
    (; N, P, Z, D,u, v, w, ζ);
                                                schedule = TimeInterval(2hours),
                                                filename = "eady_turbulence_npzd",
                                                overwrite_existing = true)
run!(simulation)

  ζ = FieldTimeSeries("eady_turbulence_npzd.jld2", "ζ")
  P = FieldTimeSeries("eady_turbulence_npzd.jld2", "P")
  N = FieldTimeSeries("eady_turbulence_npzd.jld2", "N")


times = ζ.times

xζ, yζ, zζ = nodes(ζ)
xc, yc, zc = nodes(P)

n = Observable(1)


Nz = grid.Nz

  ζₙ = @lift interior(  ζ[$n], :, :, Nz)
  Nₙ = @lift interior(N[$n], :, :, Nz)
  Pₙ = @lift interior(  P[$n], :, :, Nz)


fig = Figure(size = (1600, 1600), fontsize = 20)

lims = [(minimum(T), maximum(T)) for T in (  ζ[:, :, Nz, :],
                                           N[:, :, Nz, :] ,
                                             P[:, :, Nz, :],
                                          )]

axis_kwargs = (xlabel = "x (m)", ylabel = "y (m)", aspect = DataAspect())

ax1 = Axis(fig[1, 1]; title = "Vertical vorticity (1 / s)", axis_kwargs...)
hm1 = heatmap!(ax1, xζ, yζ, ζₙ, colormap = :balance, colorrange = lims[1])
Colorbar(fig[1, 2], hm1)

ax2 = Axis(fig[1, 3]; title = "Nutrient concentration (mmol N / m³)", axis_kwargs...)
hm2 = heatmap!(ax2, xc, yc, Nₙ, colormap = Reverse(:bamako), colorrange = lims[2])
Colorbar(fig[1, 4], hm2)

ax3 = Axis(fig[2, 1]; title = "Phytoplankton concentration (mmol N / m³)", axis_kwargs...)
hm3 = heatmap!(ax3, xc, yc, Pₙ;colorrange=(0.45,0.51), colormap=:balance)
Colorbar(fig[2, 2], hm3)



title = @lift "t = $(prettytime(times[$n]))"
Label(fig[0, :], title, fontsize = 30)


record(fig, "eady-NPZD.mp4", 2:length(times)-1, framerate=8) do i
    n[] = i
end