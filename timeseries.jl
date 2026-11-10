using Oceananigans
using CairoMakie
using Oceananigans.Units
using Statistics

# ============================================================
# 1. 参数设置
# ============================================================

filename = "Q_100_m_1_3000_res_128.jld2"

# 积分深度：上层 300 m
depth = 100.0   # m


# ============================================================
# 2. 读取三个变量的时间序列
# ============================================================

Tts = FieldTimeSeries(filename, "T")
Pts = FieldTimeSeries(filename, "P")
Nts = FieldTimeSeries(filename, "N")
w = FieldTimeSeries(filename, "w")
times = w.times

grid = Tts.grid
Nt = length(times)
Nx, Ny, Nz, Nt = size(Tts)
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
    wT′ = w′ .* T′.*1025*3990
    heat_flux[:,:,:,n]=wT′

end

# ============================================================
# 3. 获取 grid 和时间
# ============================================================



# ============================================================
# 4. 获取 z 方向的 cell faces
# ============================================================

# zf 是网格单元边界的位置
#
# 例如：
#
#       z = 0 m
#       ───────────────
#          cell 32
#       ───────────────
#          cell 31
#       ───────────────
#             ...
#       ───────────────
#          cell 1
#       ───────────────
#       z = -500 m
#
zf = znodes(grid, Face())



# ============================================================
# 5. 计算每个垂向 cell 与 [-depth, 0] 的重叠厚度
# ============================================================

Nz = length(zf) - 1

vertical_weights = zeros(Float64, Nz)

for k in 1:Nz

    # 当前 cell 的上下边界
    z_bottom = zf[k]
    z_top    = zf[k+1]

    # 当前 cell 与 [-depth, 0] 的重叠部分
    overlap_bottom = max(z_bottom, -depth)
    overlap_top    = min(z_top, 0.0)

    vertical_weights[k] =
        max(0.0, overlap_top - overlap_bottom)

end


# ============================================================
# 6. 检查实际积分深度
# ============================================================

actual_depth = sum(vertical_weights)


# ============================================================
# 7. 计算上层 300 m 的体积
# ============================================================


horizontal_area = grid.Lx * grid.Ly

volume_300m = horizontal_area * actual_depth




# ============================================================
# 8. 创建保存结果的数组
# ============================================================

nt = length(times)

T_mean = zeros(Float64, nt)
P_mean = zeros(Float64, nt)
N_mean = zeros(Float64, nt)
VHF_mean = zeros(Float64, nt)


# ============================================================
# 9. 对每一个时间帧进行积分
# ============================================================

for n in 1:nt

    println("Processing frame ", n, " / ", nt)

    # --------------------------------------------------------
    # 读取当前时间的三维场
    # --------------------------------------------------------

    Tn = Tts[n]
    Pn = Pts[n]
    Nn = Nts[n]


    # --------------------------------------------------------
    # 读取内部数组
    #
    # size 应该是：
    #
    # 128 × 128 × 32
    #
    # --------------------------------------------------------

    T_data = interior(Tn)
    P_data = interior(Pn)
    N_data = interior(Nn)


    # --------------------------------------------------------
    # 对每个 z 层进行水平平均
    # 然后按照该层在上层300m中的实际厚度加权
    # --------------------------------------------------------

    T_integral = 0.0
    P_integral = 0.0
    N_integral = 0.0


    for k in 1:Nz

        # 该层参与积分的厚度
        Δz_k = vertical_weights[k]

        # 如果这一层完全位于 300 m 以下，则跳过
        if Δz_k == 0.0
            continue
        end


        # 当前 z 层的水平平均值
        T_horizontal_mean =
            sum(@view T_data[:, :, k]) / (size(T_data, 1) * size(T_data, 2))

        P_horizontal_mean =
            sum(@view P_data[:, :, k]) / (size(P_data, 1) * size(P_data, 2))

        N_horizontal_mean =
            sum(@view N_data[:, :, k]) / (size(N_data, 1) * size(N_data, 2))


        # ----------------------------------------------------
        # 乘以该层厚度
        #
        # 因为水平面积处处相同，
        # 最后的水平面积会在分子分母中抵消。
        # ----------------------------------------------------

        T_integral += T_horizontal_mean * Δz_k

        P_integral += P_horizontal_mean * Δz_k

        N_integral += N_horizontal_mean * Δz_k

    end


    # --------------------------------------------------------
    # 除以总积分深度
    #
    # 得到上层300 m的体积平均值
    # --------------------------------------------------------

    T_mean[n] = T_integral / actual_depth

    P_mean[n] = P_integral / actual_depth

    N_mean[n] = N_integral / actual_depth

end

z = znodes(grid, Center())
z_indices = findall(z .>= -300)
for n in 1:size(heat_flux, 4)

    VHF_mean[n] =
        mean(@view heat_flux[:, :, z_indices, n])

end
# ============================================================
# 10. 转换时间单位
# ============================================================

time_hours = times ./ 3600.0

time_days = times ./ 86400.0


# ============================================================
# 11. 保存结果到 CSV 文件
# ============================================================

open("upper100m_average.csv", "w") do io

    println(io, "time_seconds,time_hours,time_days,T_mean,P_mean,N_mean")

    for n in 1:nt

        println(
            io,
            times[n], ",",
            time_hours[n], ",",
            time_days[n], ",",
            T_mean[n], ",",
            P_mean[n], ",",
            N_mean[n], ",",
            VHF_mean[n]
        )

    end

end


println()
println("Results saved to upper300m_average.csv")


# ============================================================
# 12. 画温度
# ============================================================

fig_T = Figure(size = (900, 500))

ax_T = Axis(
    fig_T[1, 1],
    xlabel = "Time (days)",
    ylabel = "Volume-averaged T",
    title = "Upper 100 m Mean Temperature"
)

lines!(
    ax_T,
    time_days,
    T_mean,
    linewidth = 2
)

save("T_upper100m_Q_0_m_1_3000_res_128.png", fig_T)


# ============================================================
# 13. 画浮游植物
# ============================================================

fig_P = Figure(size = (900, 500))

ax_P = Axis(
    fig_P[1, 1],
    xlabel = "Time (days)",
    ylabel = "Volume-averaged P",
    title = "Upper 100 m Mean Phytoplankton"
)

lines!(
    ax_P,
    time_days,
    P_mean,
    linewidth = 2

)

save("P_upper100m_Q_0_m_1_3000_res_128.png", fig_P)


# ============================================================
# 14. 画营养盐
# ============================================================

fig_N = Figure(size = (900, 500))

ax_N = Axis(
    fig_N[1, 1],
    xlabel = "Time (days)",
    ylabel = "Volume-averaged N",
    title = "Upper 100 m Mean Nutrient"
)

lines!(
    ax_N,
    time_days,
    N_mean,
    linewidth = 2
)

save("N_upper100m_Q_0_m_1_3000_res_128.png", fig_N)

# ============================================================
# 14. 画VHF
# ============================================================

fig_VHF = Figure(size = (900, 500))

ax_VHF = Axis(
    fig_VHF[1, 1],
    xlabel = "Time (days)",
    ylabel = "Volume-averaged VHF",
    title = "Upper 100 m Mean VHF"
)

lines!(
    ax_VHF,
    time_days,
    VHF_mean,
    linewidth = 2
)

save("VHF_upper100m_Q_100_m_0_res_128.png", fig_VHF)
# ============================================================
# 15. 显示图
# ============================================================


