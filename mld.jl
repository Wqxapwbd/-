
run!(simulation)
ζ = FieldTimeSeries("bottom_cool.jld2", "ζ")
N = FieldTimeSeries("bottom_cool.jld2", "N")
P = FieldTimeSeries("bottom_cool.jld2", "P")
Z = FieldTimeSeries("bottom_cool.jld2", "Z")
D = FieldTimeSeries("bottom_cool.jld2", "D")
S = FieldTimeSeries("bottom_cool.jld2", "S")
times = ζ.times

xζ, yζ, zζ = nodes(ζ)
xc, yc, zc = nodes(P)
n = Observable(1)

Nz = grid.Nz

  ζₙ = @lift interior(  ζ[$n], :, :, Nz)
  Pₙ = @lift interior(  P[$n], :, :, Nz)
  Nₙ = @lift interior(N[$n], :, :, Nz)
  Dₙ = @lift interior(D[$n], :, :, Nz)
fig = Figure(size = (1600, 1600), fontsize = 20)

lims = [(minimum(T), maximum(T)) for T in (  ζ[:, :, Nz, :],
                                           N[:, :, Nz, :] ,
                                             P[:, :, Nz, :],
                                           D[:, :, Nz, :])]

axis_kwargs = (xlabel = "x (m)", ylabel = "y (m)", aspect = DataAspect())

ax1 = Axis(fig[1, 1]; title = "Vertical vorticity (1 / s)", axis_kwargs...)
hm1 = heatmap!(ax1, xζ, yζ, ζₙ, colormap = :balance, colorrange = (-0.001,0.001))
Colorbar(fig[1, 2], hm1)

ax2 = Axis(fig[1, 3]; title = "Nutrient concentration (mmol N / m³)", axis_kwargs...)
hm2 = heatmap!(ax2, xc, yc, Nₙ, colormap = Reverse(:bamako), colorrange = lims[2])
Colorbar(fig[1, 4], hm2)

ax3 = Axis(fig[2, 1]; title = "Phytoplankton concentration (mmol N / m³)", axis_kwargs...)
hm3 = heatmap!(ax3, xc, yc, Pₙ, colormap = Reverse(:batlow), colorrange = lims[3])
Colorbar(fig[2, 2], hm3)

ax4 = Axis(fig[2, 3]; title = "Detritus concentration", axis_kwargs...)
hm4 = heatmap!(ax4, xc, yc, Dₙ, colormap = Reverse(:devon), colorrange = lims[4])
Colorbar(fig[2, 4], hm4)

title = @lift "t = $(prettytime(times[$n]))"
Label(fig[0, :], title, fontsize = 30)
record(fig, "surface cooling Q=-10 ver2.mp4", 1:length(times), framerate = 12) do i
    n[] = i
end

#读取MLD
using Oceananigans
using JLD2
using Statistics
using Printf
using CairoMakie

Tts = FieldTimeSeries(
    "surfacecooling-1.jld2",
    "T"
)
times = Tts.times
grid = Tts.grid
zc = znodes(grid, Center())

function temperature_to_density(T)

    ρ0 = 1025.0

    α = 2e-4   # 海水热膨胀系数

    T0 = 5.0   # 你的初始温度

    ρ = ρ0 .* (1 .- α .* (T .- T0))

    return ρ

end

function density_profile(T)

    ρ = temperature_to_density(T)

    ρprofile =
        dropdims(
            mean(ρ,dims=(1,2)),
            dims=(1,2)
        )

    return ρprofile

end

function calculate_MLD(
    ρprofile,
    zc;
    Δρ=0.03
)


    ρsurface = ρprofile[1]


    Δ =
        abs.(ρprofile .- ρsurface)


    iz =
        findfirst(Δ .>= Δρ)


    if isnothing(iz)

        return abs(zc[end])

    else

        return abs(zc[iz])

    end

end

MLD = zeros(length(times))


for n in eachindex(times)

    println(n)

    T =
      interior(Tts[n])


    ρprofile =
      density_profile(T)


    MLD[n] =
      calculate_MLD(
          ρprofile,
          zc
      )

end

#比较不同surface cooling的MLD
MLD1 = calculate_MLD_timeseries(
    "surfacecooling-1.jld2"
)

MLD10 = calculate_MLD_timeseries(
    "surfacecooling-10.jld2"
)

MLD100 = calculate_MLD_timeseries(
    "surfacecooling-100.jld2"
)

lines!(ax,
       times/86400,
       MLD,
       linewidth=3)


display(fig)

lines!(ax,
       times,
       MLD1,
       label="1 W/m²")

lines!(ax,
       times,
       MLD10,
       label="10 W/m²")

lines!(ax,
       times,
       MLD100,
       label="100 W/m²")


axislegend(ax)


fig = Figure(size=(500,700))


ax = Axis(
    fig[1,1],
    xlabel = "Density (kg m⁻³)",
    ylabel = "Depth (m)",
    title = "Density profile at (64,64)"
)


lines!(
    ax,
    ρ_profile,
    z,
    linewidth=3
)


ax.yreversed = true


fig