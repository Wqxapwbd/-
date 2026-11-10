P_ed=FieldTimeSeries("Ri_1_res_128_Q_300.jld2","P")
P_non=FieldTimeSeries("Ri_1_res_16_Q_300.jld2","P")
times = P_non.times
days=times/86400
Nz=32

Nx, Ny, Nz, Nt = size(P_ed)



A = interior(P_ed[1])
B = interior(P_non[1])

fig = Figure(size=(800,600))
ax = Axis(fig[1,1],
          xlabel="x",
          ylabel="y",
          title="Phytoplankton Concentration Difference")


n = Observable(1)

Pdiff = @lift begin

    P1_now = interior(P_ed_mean[$n], :, :, Nz)
    P2_now = interior(P_non_mean[$n], :, :, Nz)

    P1_now .- P2_now

end

hm = heatmap!(
    ax,
    Pdiff,
    colormap=:balance,
    colorrange=(-0.05,0.05)
)
Colorbar(
    fig[1,2],
    hm
)

title = @lift "t = $(prettytime(times[$n]))"
Label(fig[0, :], title, fontsize = 30)
record(fig,"Q_100 Pdiff.mp4",1:length(P_non.times)) do i

    n[] = i

end



#不同分辨率取体积平均的对比
P_ed_mean = zeros( Nt);
P_non_mean = zeros( Nt);

for n in 1:Nt

    # 三维数组 (Nx, Ny, Nz)
    Ped_n = interior(P_ed[n])
    Pnon_n = interior(P_non[n])


    # 水平平均
    P_ed_mean[n] = mean(Ped_n)
    P_non_mean[n] = mean(Pnon_n)

end

P_dif=P_ed_mean .- P_non_mean
P_rel=P_dif./P_non_mean

fig = Figure(size=(800,600))
ax = Axis(fig[1,1],
          xlabel="days",
          ylabel="P_dif/P_non",
          title="Phytoplankton Concentration Difference (Q=300)")

lines!(
    ax,
    days,
    P_rel,
    linewidth = 2
)

w_ed = FieldTimeSeries("Ri_1_res_128_Q_300.jld2", "w")
w_non = FieldTimeSeries("Ri_1_res_16_Q_300.jld2", "w")
Tts_ed = FieldTimeSeries("Ri_1_res_128_Q_300.jld2", "T")
Tts_non = FieldTimeSeries("Ri_1_res_16_Q_300.jld2", "T")
Nt = length(Tts_ed.times)

times = w_non.times
days=times/86400

heat_flux_ed = zeros(Float64, 32, Nt);
heat_flux_non = zeros(Float64, 32, Nt);

for n in 1:Nt

    # 三维数组 (Nx, Ny, Nz)
    Tn_ed = interior(Tts_ed[n])
    wn_ed = interior(w_ed[n])
    w_center_ed = 0.5 .* (wn_ed[:, :, 1:end-1] .+ wn_ed[:, :, 2:end])

    Tn_non = interior(Tts_non[n])
    wn_non = interior(w_non[n])
    w_center_non = 0.5 .* (wn_non[:, :, 1:end-1] .+ wn_non[:, :, 2:end])

    # 水平平均
    T̄_ed = mean(Tn_ed, dims=(1, 2))
    w̄_ed = mean(w_center_ed, dims=(1, 2))
     T̄_non = mean(Tn_non, dims=(1, 2))
    w̄_non = mean(w_center_non, dims=(1, 2))

    # 扰动
    T_ed′ = Tn_ed .- T̄_ed
    w_ed′ = w_center_ed .- w̄_ed

    T_non′ = Tn_non .- T̄_non
    w_non′ = w_center_non .- w̄_non

    # 扰动热通量
    wT_ed = @.w_ed′ .* T_ed′.*1025*3990
    heat_flux_ed[:, n] = dropdims(mean(wT_ed, dims=(1, 2)), dims=(1, 2))
    wT_non = @.w_non′ .* T_non′.*1025*3990
    heat_flux_non[:, n] = dropdims(mean(wT_non, dims=(1, 2)), dims=(1, 2))

end

vhf_dif=heat_flux_ed[30, :]-heat_flux_non[30, :]
fig = Figure(size=(800,600))
ax = Axis(fig[1,1],
          xlabel="days",
          ylabel="vhf_dif",
          title="VHF difference (Q=300)")

lines!(
    ax,
    days,
    vhf_dif,
    linewidth = 2
)

