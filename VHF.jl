using Statistics
using OceanBioME, Oceananigans, Printf, CairoMakie
using Oceananigans.Units

using Random
Random.seed!(11)
filename = "Q_0_m_0_res_128.jld2"

Tts = FieldTimeSeries(filename, "T")
w = FieldTimeSeries(filename, "w")
ζ = FieldTimeSeries(filename, "ζ")
P = FieldTimeSeries(filename, "P")
times = ζ.times
grid= w.grid

Nt = length(times)
Nz = grid.Nz
Nx=128
Ny=128
heat_flux = zeros(Float64, Nx,Ny,Nz, Nt);

for n in 1:Nt

    # 三维数组 (Nx, Ny, Nz)
    Tn = interior(Tts[n])
    wn = interior(w[n])
    w_center = 0.5 .* (wn[:, :, 1:end-1] .+ wn[:, :, 2:end])

    # 水平平均
    T̄ = mean(Tn, dims=(1, 2))
    w̄ = mean(w_center, dims=(1, 2))

    # 扰动
    T′ = Tn .- T̄
    w′ = w_center .- w̄


    # 扰动热通量
    wT′ = w′ .* T′
    heat_flux[:,:,:,n]=wT′

end

idx=64
idy=64

xc, yc, zc = nodes(P)

using CairoMakie
fig = Figure(size=(700,500))
ax = Axis(fig[1,1],
          xlabel = "w'T' (K m s⁻¹)",
          ylabel = "z (m)",
          title = "Vertical turbulent heat flux")

lines!(ax,
       heat_flux[idx,idy,:, end],
       zc)

fig

#画整个水平面的图
fig = Figure(size=(700,600))

ax = Axis(
    fig[1,1],
    xlabel="x",
    ylabel="y",
    title="VHF"
)


hm = heatmap!(
    ax,
    heat_flux[:,:,32,14],colormap = :balance, colorrange = (-0.0001,0.0001)
)
Colorbar(fig[1, 2], hm)


record(fig, "VHF_m_8e-5_z32.mp4", 1:size(heat_flux,4), framerate=10) do n

    hm[3][] = heat_flux[:,:,30,n]
    t_days = Tts.times[n] / 86400

    ax.title = "Time = $(round(t_days, digits=2)) days"

end

fig
save("VHF_with_eddy_30.png",fig)


title = @lift "t = $(prettytime(times[$n]))"
Label(fig[0, :], title, fontsize = 30)
record(fig,"VHF_with_eddy.mp4",242:length(Tts.times)) do i

    n[] = i

end

#含时间的图
fig = Figure(size=(900,600))

ax = Axis(fig[1,1],
          xlabel="Time",
          ylabel="Depth")

hm = heatmap!(ax,
              times,
              zc,
              heat_flux,
              colormap=:balance)

Colorbar(fig[1,2], hm)
fig

n = Observable(1)

fluxₙ = @lift heat_flux[:, $n]

fig = Figure(size=(800,1200))

ax = Axis(fig[1,1],
          xlabel="w′T′",
          ylabel="z",
          aspect=0.5)

rowsize!(fig.layout, 1, Relative(0.8))
xlims!(ax,-0.00001,0.00001)
lines!(ax, fluxₙ, zc)

title = @lift "t = $(times[$n])"

Label(fig[0,1], title)


record(fig,
       "VHF.mp4",
       1:length(times),
       framerate=12) do i

    n[] = i

end

using Oceananigans
using Statistics

Tts = FieldTimeSeries(
    "eady_npzd.jld2",
    "T"
)

wts = FieldTimeSeries(
    "eady_npzd.jld2",
    "w"
)


Nt = length(Tts.times)

Nx, Ny, Nz = size(interior(Tts[20]))


# 保存 x-z切面
VHF_xz = zeros(Nx, Nz, Nt)

for n in 1:Nt


    # 温度
    Tn = interior(Tts[n])


    # 垂向速度
    wn = interior(wts[n])


    # w插值到中心
    w_center =
        0.5 .* (
        wn[:, :, 1:end-1]
        .+
        wn[:, :, 2:end]
        )


    # 选择一个y切面
    T_xz =
        Tn[:, j0, :]


    w_xz =
        w_center[:, j0, :]


    # 计算水平平均背景
    Tbar =
        mean(Tn, dims=(1,2))


    wbar =
        mean(w_center, dims=(1,2))


    # 转成 Nz
    Tbar =
        vec(Tbar)


    wbar =
        vec(wbar)


    # 扰动
    Tprime =
        T_xz .-
        reshape(Tbar,1,Nz)


    wprime =
        w_xz .-
        reshape(wbar,1,Nz)


    # VHF
    VHF_xz[:,:,n] .=
        wprime .* Tprime


end
