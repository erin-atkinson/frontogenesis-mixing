function energytimeseries(run_ids)
    
    fig = Figure(; 
        fontsize = 18,
        size = (800, 800),
    )

    ax_ke = Axis(fig[1, 1]; xlabel=t_label, ylabel=L"\text{Kinetic energy} / \rho L_D^2 f^2")
    ax_pe = Axis(fig[1, 2]; xlabel=t_label, ylabel=L"\text{Potential energy} / \rho L_D^2 f^2")
    hidexdecorations(ax_ke; ticks=false, grid=false)
    hidexdecorations(ax_pe; ticks=false, grid=false)

    ax_ζ = Axis(fig[2, 1]; xlabel=t_label, ylabel=L"\zeta_\pm / f")
    ax_δ = Axis(fig[2, 2]; xlabel=t_label, ylabel=L"\delta_\pm / f")

    colors = Makie.wong_colors()
    
    N = length(run_ids)
    lns = map(1:N) do i
        run_id = run_ids[i]
        color = colors[i]

        ENERGY = joinpath(scratchpath, run_id, "ENERGY.jld2")
        STRAINTENSOR = joinpath(scratchpath, run_id, "STRAINTENSOR.jld2")

        ke = interior(FieldTimeSeries(ENERGY, "ke"), 1, 1, 1, :)
        pe = interior(FieldTimeSeries(ENERGY, "pe"), 1, 1, 1, :)

        ζ = FieldTimeSeries(STRAINTENSOR, "ζ")
        δ = FieldTimeSeries(STRAINTENSOR, "δ")

        ζ_pos = similar(times)
        ζ_neg = similar(times)

        δ_pos = similar(times)
        δ_neg = similar(times)

        for i in 1:length(times)
            ζi = interior(ζ, :, :, 1, i)
            δi = interior(δ, :, :, 1, i)

            (ζ₊, n₊, ζ₋, n₋) = mapreduce(a->accsigns(a, 0.01 * sp.f), .+, ζi)
            ζ_pos[i] = ζ₊ / max(n₊, 1)
            ζ_neg[i] = ζ₋ / max(n₋, 1)

            (δ₊, n₊, δ₋, n₋) = mapreduce(a->accsigns(a, 0.01 * sp.f), .+, δi)
            δ_pos[i] = δ₊ / max(n₊, 1)
            δ_neg[i] = δ₋ / max(n₋, 1)
        end

        lines!(ax_ke, times ./ t_unit, ke ./ (sp.L * sp.f)^2; color)
        lines!(ax_pe, times ./ t_unit, pe ./ (sp.L * sp.f)^2; color)

        lines!(ax_ζ, times ./ t_unit, ζ_pos ./ sp.f; color)
        lines!(ax_ζ, times ./ t_unit, ζ_neg ./ sp.f; color)

        lines!(ax_δ, times ./ t_unit, δ_pos ./ sp.f; color)
        lines!(ax_δ, times ./ t_unit, δ_neg ./ sp.f; color)
    end

    Legend(fig[1:2, 3], lns, run_ids)
    return fig
end

function accsigns(ζ, threshold=zero(ζ))
    ζ > threshold && return (ζ, 1, zero(ζ), 0)
    ζ < -threshold && return (zero(ζ), 0, ζ, 1)
    return (zero(ζ), 0, zero(ζ), 0)
end
