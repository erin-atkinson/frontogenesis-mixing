function frontogenesisspectra(run_id)
    FRONTOGENESIS = joinpath(scratchpath, run_id, "FRONTOGENESIS.jld2")
    iterations, times = iterations_times(FRONTOGENESIS)
    iteration = iterations[end]
    sp = simulation_parameters(FRONTOGENESIS)

    mixing = get_field(FRONTOGENESIS, "mixing_x", iteration) .+ get_field(FRONTOGENESIS, "mixing_y", iteration)
    advection = get_field(FRONTOGENESIS, "advection_x", iteration) .+ get_field(FRONTOGENESIS, "advection_y", iteration)
    Fh = get_field(FRONTOGENESIS, "Fh", iteration)
    b_mld = get_field(FRONTOGENESIS, "b_mld", iteration)

    ks, mixing_fft = modedecomposition(interior(mixing, :, :, 1), interior(b_mld, :, :, 1), sp.Lx/sp.Nx)
    ks, advection_fft = modedecomposition(interior(advection, :, :, 1), interior(b_mld, :, :, 1), sp.Lx/sp.Nx)
    ks, Fh_fft = modedecomposition(interior(Fh, :, :, 1), interior(b_mld, :, :, 1), sp.Lx/sp.Nx)

    fig = Figure()
    ax = Axis(fig[1, 1];
        xlabel = L"k",
        ylabel = L"\text{Contribution}"
    )

    lines!(ax, ks, mixing_fft)
    lines!(ax, ks, advection_fft)
    lines!(ax, ks, Fh_fft)

    fig
end
