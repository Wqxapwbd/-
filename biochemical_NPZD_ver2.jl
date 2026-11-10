using OceanBioME, Oceananigans, Printf, CairoMakie
using Oceananigans.Units

using Random
Random.seed!(11)
grid = RectilinearGrid(size = (128,128, 32),
                       extent = (10kilometer, 10kilometer, 500meters))
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
biogeochemistry = NPZD(grid;
                          scale_negatives = true)

model = NonhydrostaticModel(grid;
                            biogeochemistry,
                            advection = WENO(),
                            coriolis,
                            tracers = (:N, :P, :Z,:T,:S ,:D),
                            buoyancy,
                            background_fields = (T = T_field, v = V_field),
                            closure = vertical_diffusivity)
Ξ(z) = randn() * z / grid.Lz * (z / grid.Lz + 1)

Ũ = 1e-3
uᵢ(x, y, z) = Ũ * Ξ(z)
vᵢ(x, y, z) = Ũ * Ξ(z)

set!(model, u=uᵢ, v=vᵢ, P = 0.03, Z = 0.03, N = 4.0, S = 35, T = 20)
Δx = minimum_xspacing(grid, Center(), Center(), Center())
Δy = minimum_yspacing(grid, Center(), Center(), Center())
Δz = minimum_zspacing(grid, Center(), Center(), Center())

Δt₀ = 0.5 * min(Δx, Δy, Δz) / V(0, 0, 0, 0, background_state_parameters)

simulation = Simulation(model, Δt = Δt₀, stop_time = 10days)
wizard = TimeStepWizard(cfl = 0.5, diffusive_cfl = 0.5, max_Δt = 20minutes)
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
simulation.output_writers[:fields] = JLD2Writer(model, merge(model.tracers, (; u, v, w, ζ));
                                                schedule = TimeInterval(2hours),
                                                filename = "non_eddy_2",
                                                overwrite_existing = true)
run!(simulation)

filename="Q_100_m_1_3000_res_128.jld2"

ζ = FieldTimeSeries(filename, "ζ")
N = FieldTimeSeries(filename, "N")
P = FieldTimeSeries(filename, "P")
Z = FieldTimeSeries(filename, "Z")
D = FieldTimeSeries(filename, "D")
Tt= FieldTimeSeries(filename, "D")

times = ζ.times
grid = N.grid

xζ, yζ, zζ = nodes(ζ)
xc, yc, zc = nodes(P)
n = Observable(1)

Nz = grid.Nz

  ζₙ = @lift interior(  ζ[$n], :, :, 32)
  Pₙ = @lift interior(  P[$n], :, :, 32)
  Nₙ = @lift interior(N[$n], :, :, 32)
  Dₙ = @lift interior(D[$n], :, :, 32)
  Tn=@lift interior(Tt[$n], :, :, 32)
fig = Figure(size = (1600, 1600), fontsize = 20)

lims = [(minimum(T), maximum(T)) for T in (  ζ[:, :, Nz, :],
                                           N[:, :, Nz, :] ,
                                             P[:, :, Nz, :],
                                           D[:, :, Nz, :])]


axis_kwargs = (xlabel = "x (m)", ylabel = "y (m)", aspect = DataAspect())

ax1 = Axis(fig[1, 1]; title = "Vertical vorticity (1 / s)", axis_kwargs...)
hm1 = heatmap!(ax1, xζ, yζ, ζₙ, colormap = :balance, colorrange = (-0.002,0.002))
Colorbar(fig[1, 2], hm1)

ax2 = Axis(fig[1, 3]; title = "Nutrient concentration (mmol N / m³)", axis_kwargs...)
hm2 = heatmap!(ax2, xc, yc, Nₙ, colormap = Reverse(:bamako), colorrange = lims[2])
Colorbar(fig[1, 4], hm2)

ax3 = Axis(fig[2, 1]; title = "Phytoplankton concentration relative difference(mmol N / m³)", axis_kwargs...)
hm3 = heatmap!(ax3, xc, yc, Pₙ, colormap = Reverse(:batlow), colorrange =lims[3])
Colorbar(fig[2, 2], hm3)

ax4 = Axis(fig[2, 3]; title = "Detritus concentration", axis_kwargs...)
hm4 = heatmap!(ax4, xc, yc, Dₙ, colormap = Reverse(:devon), colorrange = lims[4])
Colorbar(fig[2, 4], hm4)

#新增温度图
lims = [(minimum(T), maximum(T)) for T in (Tt[:, :, Nz, :])]
ax4 = Axis(fig[2, 3]; title = "Temp", axis_kwargs...)
hm4 = heatmap!(ax4, xc, yc, Tn, colormap = Reverse(:devon), colorrange = lims[1])
Colorbar(fig[2, 4], hm4)

title = @lift "t = $(prettytime(times[$n]))"
Label(fig[0, :], title, fontsize = 30)
record(fig, "Q_100_m_1_3000_res_128.mp4", 1:length(times), framerate = 12) do i
    n[] = i
end
