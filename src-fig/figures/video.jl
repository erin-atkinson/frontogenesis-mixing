function gradient_video(run_id, frames, filename;
    fig_kw = NamedTuple(),
    record_kw = NamedTuple(),
    )
    
    foldername = joinpath(scratchpath, run_id)
    
    SURFACE = joinpath(foldername, "SURFACE.jld2")

    sp = simulation_parameters(SURFACE)
    iterations, times = iterations_times(SURFACE)
    
    n = Observable(frames[1])
    t = @lift interp_time($n, times)
    
    fts_u = FieldTimeSeries(SURFACE, "u"; backend=OnDisk())
    fts_v = FieldTimeSeries(SURFACE, "v"; backend=OnDisk())
    #fts_b = FieldTimeSeries(SURFACE, "b"; backend=OnDisk())

    temp_u = similar(fts_u[1])
    temp_v = similar(fts_v[1])
    
    vorticity = Field(∂x(temp_v) - ∂y(temp_u))
    divergence = Field(∂x(temp_u) + ∂y(temp_v))
    strain = Field((∂x(temp_u) - ∂y(temp_v))^2 + (∂x(temp_v) + ∂y(temp_u))^2)

    ζ = Observable(nov(vorticity[:, :, 1]) ./ sp.f)
    δ = Observable(nov(divergence[:, :, 1]) ./ sp.f)
    σ = Observable(nov(strain[:, :, 1]) ./ sp.f^2)
    
    on(t) do t
        set!(temp_u, fts_u[Time(t)])
        set!(temp_v, fts_v[Time(t)])
        
        compute!(vorticity)
        compute!(divergence)
        compute!(strain)
        
        ζ[] = nov(vorticity[:, :, 1]) ./ sp.f
        δ[] = nov(divergence[:, :, 1]) ./ sp.f
        σ[] = nov(strain[:, :, 1]) ./ sp.f^2
    end
    
    title = @lift begin
        t_day = @sprintf "%.0f" ($t / (24*3600))
        ft_str = @sprintf "%.1f" sp.f * $t
        fT_mix_str = @sprintf "%.1f" sp.f * sp.T_mix
        L"\text{Submesoscale} \quad fT_\text{mix} = %$fT_mix_str \quad t = %$t_day \, \text{day}\quad ft = %$ft_str"
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
        heatmap!(ax_ζ, xs, ys, data; colormap=:curl, colorrange=(-1, 1))
    end

    ht_δ = begin
        xs = nov(xnodes(divergence; with_halos=true)) ./ sp.L
        ys = nov(ynodes(divergence; with_halos=true)) ./ sp.L
        data = δ
        heatmap!(ax_δ, xs, ys, data; colormap=:balance, colorrange=(-1, 1))
    end

    ht_σ = begin
        xs = nov(xnodes(strain; with_halos=true)) ./ sp.L
        ys = nov(ynodes(strain; with_halos=true)) ./ sp.L
        data = σ
        heatmap!(ax_σ, xs, ys, data; colormap=:balance, colorrange=(-1, 1))
    end

    Colorbar(fig[3, 1], ht_ζ; label=L"\zeta / f", vertical=false, flipaxis=false)
    Colorbar(fig[3, 2], ht_δ; label=L"\delta / f", vertical=false, flipaxis=false)
    Colorbar(fig[3, 3], ht_σ; label=L"\sigma / f", vertical=false, flipaxis=false)

    colgap!(fig.layout, 40)
    prettyrecord(n, fig, filename, frames; record_kw...)

    return fig
end
