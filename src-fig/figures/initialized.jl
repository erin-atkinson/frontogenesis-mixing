function initialized_BI(run_id, frames, filename;
    fig_kw = NamedTuple(),
    record_kw = NamedTuple(),
    )
    
    foldername = joinpath(initpath, run_id)
    
    INS = joinpath(foldername, "AVG.jld2")

    sp = simulation_parameters(INS)
    iterations, times = iterations_times(INS)
    
    n = Observable(frames[1])
    t = @lift interp_time($n, times)
    
    fts_u = FieldTimeSeries(INS, "u"; backend=OnDisk())
    fts_v = FieldTimeSeries(INS, "v"; backend=OnDisk())
    #fts_b = FieldTimeSeries(INS, "b"; backend=OnDisk())

    temp_u = similar(fts_u[1])
    temp_v = similar(fts_v[1])
    
    vorticity = Field(∂x(temp_v) - ∂y(temp_u))
    divergence = Field(∂x(temp_u) + ∂y(temp_v))
    strain = Field(sqrt((∂x(temp_u) - ∂y(temp_v))^2 + (∂y(temp_u) + ∂x(temp_v))^2))

    ζ = Observable(nov(vorticity[:, :, sp.Nz]) ./ sp.f)
    δ = Observable(nov(divergence[:, :, sp.Nz]) ./ sp.f)
    σ = Observable(nov(strain[:, :, sp.Nz]) ./ sp.f)
    
    on(t) do t
        set!(temp_u, fts_u[Time(t)])
        set!(temp_v, fts_v[Time(t)])
        
        compute!(vorticity)
        compute!(divergence)
        compute!(strain)
        
        ζ[] = nov(vorticity[:, :, sp.Nz]) ./ sp.f
        δ[] = nov(divergence[:, :, sp.Nz]) ./ sp.f
        σ[] = abs.(nov(strain[:, :, sp.Nz])) ./ sp.f
    end
    
    title = @lift begin
        t_day = @sprintf "%.0f" ($t / (24*3600))
        σt_str = @sprintf "%.1f" sp.σ * $t
        L"\text{Mesoscale initialisation} \quad t = %$t_day \, \text{day}\quad \sigma_\text{E}t = %$σt_str"
    end
    
    #u = @lift nov(fts_u[Time($t)][:, 1, :])
    #v = @lift nov(fts_v[Time($t)][:, 1, :])
    #b = @lift nov(fts_b[Time($t)][:, :, sp.Nz])

    fig = Figure(; 
        size=(800, 400),
        fontsize = 18,
        fig_kw...
    )
    Label(fig[1, 1:3], title)
    
    ax_kw = (;
        xlabel = L"x / L_D",
        ylabel = L"y / L_D",
        limits = (-sp.Lx / 2sp.L, sp.Lx / 2sp.L, -sp.Ly / 2sp.L, sp.Ly / 2sp.L),
    )

    ax_ζ = Axis(fig[2, 1]; ax_kw...)
    ax_δ = Axis(fig[2, 2]; ax_kw...)
    ax_σ = Axis(fig[2, 3]; ax_kw...)

    hideydecorations!(ax_δ; ticks=false)
    hideydecorations!(ax_σ; ticks=false)

    ht_ζ = begin
        xs = nov(xnodes(vorticity; with_halos=true)) ./ sp.L
        ys = nov(ynodes(vorticity; with_halos=true)) ./ sp.L
        data = ζ
        heatmap!(ax_ζ, xs, ys, data; colormap=:curl, colorrange=(-0.1, 0.1))
    end

    ht_δ = begin
        xs = nov(xnodes(divergence; with_halos=true)) ./ sp.L
        ys = nov(ynodes(divergence; with_halos=true)) ./ sp.L
        data = δ
        heatmap!(ax_δ, xs, ys, data; colormap=:balance, colorrange=(-0.1, 0.1))
    end

    ht_σ = begin
        xs = nov(xnodes(strain; with_halos=true)) ./ sp.L
        ys = nov(ynodes(strain; with_halos=true)) ./ sp.L
        data = σ
        heatmap!(ax_σ, xs, ys, data; colormap=:amp, colorrange=(0, 0.1))
    end

    Colorbar(fig[3, 1], ht_ζ; label=L"\zeta / f", vertical=false, flipaxis=false)
    Colorbar(fig[3, 2], ht_δ; label=L"\delta / f", vertical=false, flipaxis=false)
    Colorbar(fig[3, 3], ht_σ; label=L"\sigma / f", vertical=false, flipaxis=false)

    colgap!(fig.layout, 40)
    prettyrecord(n, fig, filename, frames; record_kw...)

    return fig
end