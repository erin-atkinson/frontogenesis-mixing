function energytimeseries(run_ids)
    
    fig = Figure(; 
        fontsize = 18,
        size = (1000, 600),
    )

    ax_ke = Axis(fig[1, 1]; xlabel=t_label, ylabel=L"\text{Kinetic energy} / \rho V L_D^2 f^2")

    ax_ζ = Axis(fig[2, 1]; xlabel=t_label, ylabel=L"\zeta_\pm / f")
    ax_δ = Axis(fig[2, 2]; xlabel=t_label, ylabel=L"\delta_\pm / f")

    colors = Makie.wong_colors()
    
    N = length(run_ids)
    lns = map(1:N) do r
        run_id = run_ids[r]
        color = colors[r]

        ENERGY = joinpath(scratchpath, run_id, "ENERGY.jld2")
        STRAINTENSOR = joinpath(scratchpath, run_id, "STRAINTENSOR.jld2")

        ke = interior(FieldTimeSeries(ENERGY, "ke"), 1, 1, 1, :)

        ζ = FieldTimeSeries(STRAINTENSOR, "ζ"; backend=OnDisk())
        δ = FieldTimeSeries(STRAINTENSOR, "δ"; backend=OnDisk())

        ketimes = FieldTimeSeries(ENERGY, "ke").times
        straintimes = ζ.times

        ζ_pos = similar(ketimes)
        ζ_neg = similar(ketimes)

        δ_pos = similar(ketimes)
        δ_neg = similar(ketimes)
        #=
        for i in 1:length(ketimes)
            print("$r: $i\r")
            ζi = interior(ζ[Time(ketimes[i])], :, :, 1)[:]
            δi = interior(δ[Time(ketimes[i])], :, :, 1)[:]

            (ζ₊, n₊, ζ₋, n₋) = foldr((ζ, a)->accsigns(ζ, 0.01 * sp.f) .+ a, ζi)
            ζ_pos[i] = ζ₊ / max(n₊, 1)
            ζ_neg[i] = ζ₋ / max(n₋, 1)

            (δ₊, n₊, δ₋, n₋) = foldr((δ, a)->accsigns(δ, 0.01 * sp.f) .+ a, δi)
            δ_pos[i] = δ₊ / max(n₊, 1)
            δ_neg[i] = δ₋ / max(n₋, 1)
        end
        =#
        for i in 1:length(ketimes)
            print("$r: $i\r")
            ζi = interior(ζ[Time(ketimes[i])], :, :, 1)[:]
            δi = interior(δ[Time(ketimes[i])], :, :, 1)[:]

            ζ_pos[i] = maximum(ζi)
            ζ_neg[i] = minimum(ζi)

            δ_pos[i] = maximum(δi)
            δ_neg[i] = minimum(δi)
        end

        V = sp.Lx * sp.Ly * sp.Lz
        lines!(ax_ke, ketimes ./ t_unit, ke ./ (sp.L * sp.f)^2 / V; color)

        lines!(ax_ζ, ketimes ./ t_unit, ζ_pos ./ sp.f; color)
        lines!(ax_ζ, ketimes ./ t_unit, ζ_neg ./ sp.f; color)

        lines!(ax_δ, ketimes ./ t_unit, δ_pos ./ sp.f; color)
        lines!(ax_δ, ketimes ./ t_unit, δ_neg ./ sp.f; color)
    end

    Legend(fig[1, 2], lns, run_ids; tellwidth=false)
    return fig
end

function accsigns(ζ, threshold=zero(ζ))
    ζ > threshold && return (ζ, 1, zero(ζ), 0)
    ζ < -threshold && return (zero(ζ), 0, ζ, 1)
    return (zero(ζ), 0, zero(ζ), 0)
end
