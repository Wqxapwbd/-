#画浮游植物和速度（v，u,w)的切面图
ζ = FieldTimeSeries("Ri_1_res_16_Q_300.jld2", "ζ")
P = FieldTimeSeries("Ri_1_res_16_Q_300.jld2", "P")

times = ζ.times

xζ, yζ, zζ = nodes(ζ)
xc, yc, zc = nodes(P)
n = Observable(1)
grid=P.grid
Nz = grid.Nz
Ny = 128
Nx = 128

  ζₙ = @lift interior(  ζ[$n], :, :, Nz)
  Pxz = @lift interior(  P[$n], :, Ny ,:)
  Pyz = @lift interior(  P[$n], Nx ,:,  :)
fig = Figure(size = (1600, 1600), fontsize = 20)

lims = [(minimum(T), maximum(T)) for T in (  P[:, Ny,:,  :],
                                           P[Nx, :, :, :] )]


axis_kwargs_1 = (xlabel = "x (m)", ylabel = "z (m)", aspect=6)
axis_kwargs_2 = (xlabel = "y (m)", ylabel = "z (m)",aspect=6 )

ax1 = Axis(fig[1, 1]; title = "Phytoplankton concentration P(x,z) (mmol N / m³)", axis_kwargs_1...)
hm1 = heatmap!(ax1, xc, zc, Pxz, colormap = :algae, colorrange = lims[1])
Colorbar(fig[1, 2], hm1)

ax2 = Axis(fig[2, 1]; title = "Phytoplankton concentration P(y,z) (mmol N / m³)", axis_kwargs_2...)
hm2 = heatmap!(ax2, yc, zc, Pyz, colormap = :algae, colorrange = lims[2])
Colorbar(fig[2, 2], hm2)

title = @lift "t = $(prettytime(times[$n]))"
Label(fig[0, :], title, fontsize = 30)
record(fig, "bottom_friction_P_vertical.mp4", 1:length(times), framerate = 12) do i
    n[] = i
end

filename="Q_0_m_1_3000_res_128"
u = FieldTimeSeries(filename, "u")
v = FieldTimeSeries(filename, "v")
w = FieldTimeSeries(filename, "w")
xζ, yζ, zζ = nodes(w)
xc, yc, zc = nodes(w)
n = Observable(1)
grid=w.grid
Nz = grid.Nz
uyz = @lift interior(u[$n], Nx, :, :)
vxz = @lift interior(v[$n],:, Ny, :)
wxz = @lift interior(w[$n], :, Ny, :)

fig = Figure(size = (1600, 2200), fontsize = 20)

lim = [(minimum(T), maximum(T)) for T in (  u[Nx, :,:,  :],
                                           v[:,Ny, :, :] ,
                                           w[:, Ny,:,  :] )]


axis_kwargs_1 = (xlabel = "x (m)", ylabel = "z (m)", aspect=6)
axis_kwargs_2 = (xlabel = "y (m)", ylabel = "z (m)",aspect=6 )

ax1 = Axis(fig[1, 1]; title = "U(y,z) (mmol N / m³)", axis_kwargs_2...)
hm1 = heatmap!(ax1, yc, zc, uyz, colormap = :balance, colorrange = (-0.15,0.15))
Colorbar(fig[1, 2], hm1)

ax2 = Axis(fig[2, 1]; title = "V(x,z) (mmol N / m³)", axis_kwargs_1...)
hm2 = heatmap!(ax2, xc, zc, vxz, colormap = :balance, colorrange = (-0.15,0.15))
Colorbar(fig[2, 2], hm1)

ax3 = Axis(fig[3, 1]; title = "W(y,z) (mmol N / m³)", axis_kwargs_2...)
hm3 = heatmap!(ax3, yc, zc, wxz, colormap = :balance, colorrange = (-0.15,0.15))
Colorbar(fig[3, 2], hm1)

title = @lift "t = $(prettytime(times[$n]))"
Label(fig[0, :], title, fontsize = 30)
record(fig, "Q_0_m_1_3000_res_128_vertical_uvw.mp4", 1:length(times), framerate = 12) do i
    n[] = i
end


#画密度随深度变化图
filename = "test_temp_1.jld2"

Tts = FieldTimeSeries(filename, "T")
Sts = FieldTimeSeries(filename, "S")
Tt = interior(Tts[1])
Ss = interior(Sts[1])


Ny = 128
Nx = 128
ρ0 = 1025.0

α = 2e-4      # thermal expansion coefficient
β = 7.6e-4    # haline contraction coefficient
T0 = 5.0
S0 = 35.0

rho = ρ0 .* (
    1 .-
    α .* (Tt .- T0)
    .+
    β .* (Ss .- S0)
)

ix = 32
iy = 100
ρ_profile = rho[ix, iy, :]

ρ_profile = vec(mean(rho, dims=(1,2)))

ρ_profile_ex=reverse(ρ_profile.-1000)
H = 500
Nz = 32
z = range(
    -H/32/2,
    -H + H/32/2,
    length=32
)

fig = Figure(size=(500,700))


ax = Axis(
    fig[1,1],
    xlabel = "Density(kg m⁻³)",
    ylabel = "Depth (m)",
    title = "Average Density profile "
)


lines!(
    ax,
    ρ_profile_ex,
    z,
    linewidth=3
)


ax.yreversed = false
fig

#画密度的y截面图
using CairoMakie  # 如果在 Jupyter 中或需要保存图片；如需交互窗口改用 GLMakie

# 1. 假设你的数据 rho 已经存在，尺寸为 (128, 128, 32)
# rho = rand(128, 128, 32)  # 测试用随机数据
Nx, Ny, Nz = size(Tt)
y_idx = 8

# 2. 提取 y=64 的截面
# 结果 rho_slice 的尺寸为 (128, 32)，即 (x, z)
rho_slice = Tt[:, y_idx, :]

# 3. 创建图形和坐标轴
fig = Figure(size = (800, 400)) # 设置图片宽800，高400

# 创建坐标轴，设置标签和标题
ax = Axis(fig[1, 1],
    xlabel = "x",
    ylabel = "z",
    title = "Temperature at y = $y_idx"
)

# 4. 绘制热力图 (Heatmap)
# 推荐使用 colormap 设置颜色主题，比如 :viridis, :turbo, :inferno 等
hm = heatmap!(ax, rho_slice, colormap = :turbo)

# 5. 添加颜色条 (Colorbar)
Colorbar(fig[1, 2], hm, label = "T")

# 6. (可选) 调整 z 轴方向
# 矩阵绘图时有时 z 轴会反转，如果发现 z 轴上下颠倒，取消下面这行的注释：
# ax.yreversed = true

# 7. 保存图片
save("Temp_y_64_res_128_Q_200_Ri_1.png", fig)

# 显示图片（如果在 REPL 或支持的内联环境中）
fig


using Statistics
using OceanBioME, Oceananigans, Printf, CairoMakie
using Oceananigans.Units
using Random
Random.seed!(11)

P = FieldTimeSeries("Ri_1_res_16_Q_200.jld2", "P")

# P 的大小
Nx, Ny, Nz, Nt = size(P)
grid=P.grid
# 第 n 帧
n = 121

# z 坐标
z = znodes(grid, Center())

# 每一层进行水平平均
P_profile = zeros(Float64, Nz)

for k in 1:Nz
    P_profile[k] = mean(@view P[:, :, k, n])
end

# 画图
fig = Figure(size = (600, 700))

ax = Axis(
    fig[1, 1],
    xlabel = "Phytoplankton concentration",
    ylabel = "Depth (m)",
    title = "Phytoplankton Vertical Profile",
    yreversed = false
)

lines!(
    ax,
    P_profile,
    z,
    linewidth = 2
)

scatter!(
    ax,
    P_profile,
    z
)

display(fig)


Tts = FieldTimeSeries("high_Q.jld2", "T")

# P 的大小
Nx, Ny, Nz, Nt = size(Tts)

# 第 n 帧
n = 121

# z 坐标
z = znodes(grid, Center())

# 每一层进行水平平均
Tts_profile = zeros(Float64, Nz)

for k in 1:Nz
    Tts_profile[k] = mean(@view Tts[:, :, k, n])
end

# 画图
fig = Figure(size = (600, 700))

ax = Axis(
    fig[1, 1],
    xlabel = "Temperature",
    ylabel = "Depth (m)",
    title = "Temperature Vertical Profile",
    yreversed = false
)

lines!(
    ax,
    Tts_profile,
    z,
    linewidth = 2
)

scatter!(
    ax,
    Tts_profile,
    z
)

display(fig)
