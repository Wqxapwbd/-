using OceanBioME, Oceananigans, Printf, CairoMakie
using Oceananigans.Units

using Random
Random.seed!(11)

grid = RectilinearGrid(
    size = (64, 64, 32),
    x = (0, 0kilometer),
    y = (0, 30kilometer),
    z = (-500meters, 0),
    topology = (Bounded, Periodic, Bounded)
)

# ============================================================
# 1. 设置温度剖面的参数
# ============================================================

T_surface = 7.9       # 海表温度，单位：°C
T_max = 10.0           # 最大温度，单位：°C
D_max = 100.0          # 最大温度所在深度，单位：m，输入正数

# ============================================================
# 2. 定义温度随深度变化的函数
#    z 使用海洋模型的垂向坐标：
#    海表 z = 0，海面以下 z < 0
# ============================================================

function temperature_profile(
    z;
    T_surface,
    T_max,
    D_max
)
    # 检查最大温度深度
    D_max > 0 || throw(ArgumentError("D_max 必须大于 0"))

    # 将垂向坐标 z 转换为水下深度 d（正值）
    d = -z

    # 海表以上：保持海表温度
    if d < 0
        return T_surface
    end

    # 最大温度深度以下：保持最大温度
    if d >= D_max
        return T_max
    end

    # 海表至最大温度深度之间：三次平滑插值
    s = d / D_max
    smooth = 3s^2 - 2s^3

    return T_surface + (T_max - T_surface) * smooth
end

T_initial = (x, y, z) -> temperature_profile(
    z;
    T_surface = T_surface,
    T_max = T_max,
    D_max = D_max
)


# ============================================================
# 1. 设置盐度剖面的参数
# ============================================================

S_surface = 34.0      # 海表盐度
S_max = 35.0          # 最大盐度
D_Smax = 100.0        # 最大盐度所在深度，单位：m，输入正数

# ============================================================
# 2. 定义盐度随深度变化的函数
#    海表 z = 0，海面以下 z < 0
# ============================================================

function salinity_profile(
    z;
    S_surface,
    S_max,
    D_Smax
)
    # 检查最大盐度深度
    D_Smax > 0 ||
        throw(ArgumentError("D_Smax 必须大于 0"))

    # 将垂向坐标转换为水下深度
    d = -z

    # 海表以上：保持海表盐度
    if d < 0
        return S_surface
    end

    # 最大盐度深度以下：保持最大盐度
    if d >= D_Smax
        return S_max
    end

    # 海表至最大盐度深度之间：三次平滑插值
    s = d / D_Smax
    smooth = 3s^2 - 2s^3

    return S_surface + (S_max - S_surface) * smooth
end

# 水平盐度差：x=Lx 处减去 x=0 处的盐度
ΔS_x = 0.2

# 盐度水平梯度，单位为盐度单位 / m
G_S = ΔS_x / grid.Lx

S_initial = (x, y, z) -> salinity_profile(
    z;
    S_surface = S_surface,
    S_max = S_max,
    D_Smax = D_Smax
)+ G_S * (x - grid.Lx / 2)


#记得改文件名

coriolis = FPlane(f = 1e-4) # [s⁻¹]
buoyancy = SeawaterBuoyancy()

eos = buoyancy.equation_of_state

α = eos.thermal_expansion
β = eos.haline_contraction

g = buoyancy.gravitational_acceleration
f = coriolis.f

shear = -g * β * G_S / f

@show α β g f shear


V_initial(x, y, z) = shear * (z + grid.Lz)

νᵥ = κᵥ = 1e-4 # [m² s⁻¹]
vertical_diffusivity = VerticalScalarDiffusivity(ν = νᵥ, κ = κᵥ)
biogeochemistry = NPZD(grid;
                          scale_negatives = true)

#加surface cooling和bottom drag
Qcool = 200 # W/m² 向上为正
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

                            closure = vertical_diffusivity,

                            boundary_conditions = boundary_conditions
                            )
Ξ(z) = randn() * z / grid.Lz * (z / grid.Lz + 1)

Ũ = 1e-3
uᵢ(x, y, z) = Ũ * Ξ(z)
vᵢ(x, y, z) = Ũ * Ξ(z)

set!(model, u=uᵢ, v = (x, y, z) -> V_initial(x, y, z) + vᵢ(x, y, z), P = 0.03, Z = 0.03, N = 4.0, S = S_initial, T = T_initial)
Δx = minimum_xspacing(grid, Center(), Center(), Center())
Δy = minimum_yspacing(grid, Center(), Center(), Center())
Δz = minimum_zspacing(grid, Center(), Center(), Center())

V_max = maximum(
    abs(V_initial(0, 0, z) + vᵢ(0, 0, z))
    for z in znodes(grid, Center())
)

Δt₀ = 0.5 * min(Δx, Δy, Δz) / V_max

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
                                                filename = "test_temp_large_scale.jld2",
                                                overwrite_existing = true)
run!(simulation)


#检验程序
using Statistics
using CairoMakie

# ============================================================
# 1. 提取初始温度、盐度场
# ============================================================

T_data = Array(interior(model.tracers.T))
S_data = Array(interior(model.tracers.S))

# 网格中心坐标
x = xnodes(grid, Center())
z = znodes(grid, Center())

# 水下深度：海表为 0，向下为正
depth = -z

# 温度、盐度数据的维度为 (Nx, Ny, Nz)
Nx, Ny, Nz = size(T_data)

@show size(T_data)
@show size(S_data)

# ============================================================
# 2. 检验温度是否只随深度变化
# ============================================================

# 每个深度层的水平平均温度
T_profile = vec(dropdims(mean(T_data, dims=(1, 2)), dims=(1, 2)))

# 每个深度层的水平平均盐度
S_profile = vec(dropdims(mean(S_data, dims=(1, 2)), dims=(1, 2)))

# 检查温度、盐度在每个深度层的水平变化幅度
T_horizontal_range = [
    maximum(T_data[:, :, k]) - minimum(T_data[:, :, k])
    for k in 1:Nz
]

S_horizontal_range = [
    maximum(S_data[:, :, k]) - minimum(S_data[:, :, k])
    for k in 1:Nz
]

println("温度最大水平变化幅度 = ",
        maximum(T_horizontal_range), " °C")

println("盐度最大水平变化幅度 = ",
        maximum(S_horizontal_range))

# ============================================================
# 3. 用当前线性状态方程检验水平密度梯度
# ============================================================

eos = buoyancy.equation_of_state

α = eos.thermal_expansion
β = eos.haline_contraction

g = buoyancy.gravitational_acceleration
f = coriolis.f

ρref = 1025.0

# 选取 50 m 深度附近的网格层
target_depth = 50.0
k = argmin(abs.(depth .- target_depth))

# 在该深度沿 y 方向平均，得到随 x 变化的温盐剖面
T_x = vec(mean(T_data[:, :, k], dims=2))
S_x = vec(mean(S_data[:, :, k], dims=2))

# 线性状态方程对应的相对密度异常代理
# 这里只用于检查密度梯度，不代表绝对海水密度
ρ_proxy = -ρref .* (α .* T_x .- β .* S_x)

# 去除水平平均值，便于观察密度随 x 的变化
ρ_anomaly = ρ_proxy .- mean(ρ_proxy)

# 对 x 方向数据进行线性拟合，估计密度梯度
x_centered = x .- mean(x)
ρ_centered = ρ_proxy .- mean(ρ_proxy)

dρ_dx = sum(x_centered .* ρ_centered) /
        sum(x_centered .^ 2)

# 理论上，由盐度梯度产生的密度梯度
dρ_dx_expected = ρref * β * G_S

# 由实际密度梯度反推热成风切变
shear_from_density = -g * dρ_dx / (ρref * f)

println("\n--- 水平密度梯度检验 ---")
println("检验深度 = ", depth[k], " m")
println("数值密度梯度 = ", dρ_dx, " kg/m^4")
println("线性状态方程预期梯度 = ",
        dρ_dx_expected, " kg/m^4")

println("你设定的 shear = ", shear, " s^-1")
println("由密度梯度反推的 shear = ",
        shear_from_density, " s^-1")

# ============================================================
# 4. 绘图
# ============================================================

fig = Figure(size=(1100, 750))

# 图 A：温度随深度的变化
ax1 = Axis(
    fig[1, 1],
    xlabel = "Temperature (°C)",
    ylabel = "Depth (m)",
    title = "Initial temperature profile",
    yreversed = true
)

lines!(ax1, T_profile, depth, linewidth=3)
scatter!(ax1, T_profile, depth, markersize=4)

# 图 B：盐度随深度的变化
ax2 = Axis(
    fig[1, 2],
    xlabel = "Salinity",
    ylabel = "Depth (m)",
    title = "Initial salinity profile",
    yreversed = true
)

lines!(ax2, S_profile, depth, linewidth=3)
scatter!(ax2, S_profile, depth, markersize=4)

# 图 C：固定深度的水平温度分布
ax3 = Axis(
    fig[2, 1],
    xlabel = "x (km)",
    ylabel = "Temperature (°C)",
    title = "Horizontal temperature at $(round(depth[k], digits=1)) m"
)

lines!(ax3, x ./ 1000, T_x, linewidth=3)

# 图 D：固定深度的水平密度异常
ax4 = Axis(
    fig[2, 2],
    xlabel = "x (km)",
    ylabel = "Density anomaly proxy (kg/m³)",
    title = "Horizontal density gradient at $(round(depth[k], digits=1)) m"
)

lines!(ax4, x ./ 1000, ρ_anomaly, linewidth=3)
scatter!(ax4, x ./ 1000, ρ_anomaly, markersize=5)

display(fig)

save("initial_state_diagnostics.png", fig)

