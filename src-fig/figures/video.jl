using Oceananigans: fill_halo_regions!

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
    strain = Field(sqrt((∂x(temp_u) - ∂y(temp_v))^2 + (∂x(temp_v) + ∂y(temp_u))^2))

    ζ = Observable(nov(vorticity[:, :, 1]) ./ sp.f)
    δ = Observable(nov(divergence[:, :, 1]) ./ sp.f)
    σ = Observable(nov(strain[:, :, 1]) ./ sp.f)
    
    on(t) do t
        set!(temp_u, fts_u[Time(t)])
        set!(temp_v, fts_v[Time(t)])
        
        compute!(vorticity)
        compute!(divergence)
        compute!(strain)
        
        ζ[] = nov(vorticity[:, :, 1]) ./ sp.f
        δ[] = nov(divergence[:, :, 1]) ./ sp.f
        σ[] = nov(strain[:, :, 1]) ./ sp.f
    end
    
    title = @lift begin
        t_day = @sprintf "%.0f" ($t / (24*3600))
        ft_str = @sprintf "%.1f" sp.f * $t
        γ_str =  @sprintf "%.2f" (1 / sp.T_mix / sp.f)
        L"\text{Strain invariants} \quad \gamma_\text{mix} = %$γ_str \quad t = %$t_day \, \text{day}\quad ft = %$ft_str"
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
        heatmap!(ax_σ, xs, ys, data; colormap=:amp, colorrange=(0, 1))
    end

    Colorbar(fig[3, 1], ht_ζ; label=L"\zeta / f", vertical=false, flipaxis=false)
    Colorbar(fig[3, 2], ht_δ; label=L"\delta / f", vertical=false, flipaxis=false)
    Colorbar(fig[3, 3], ht_σ; label=L"\sigma / f", vertical=false, flipaxis=false)

    colgap!(fig.layout, 40)
    prettyrecord(n, fig, filename, frames; record_kw...)

    return fig
end

function k²(k, l)
    a = k^2 + l^2
    a == 0 && return Inf
    return a
end
function helmholtzdecomposition(u, v, Δx)
    N = size(u, 1)
    u = (circshift(u, (-1, 0)) .+ circshift(u, (1, 0))) ./ 2
    v = (circshift(v, (-1, 0)) .+ circshift(v, (1, 0))) ./ 2
    
    ks = fftfreq(N, Δx)
    u_fft = fft(u)
    v_fft = fft(v)
    ψ_fft = [1im * (v_fft[i, j] * ks[i] - u_fft[i, j] * ks[j]) / k²(ks[i], ks[j]) for i in 1:N, j in 1:N]
    φ_fft = [1im * (v_fft[i, j] * ks[j]  + u_fft[i, j] * ks[i]) / k²(ks[i], ks[j]) for i in 1:N, j in 1:N]
    ψ = real.(ifft(ψ_fft))
    φ = real.(ifft(φ_fft))
    return ψ, φ
end

function potential_video(run_id, frames, filename;
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
    fts_b = FieldTimeSeries(SURFACE, "b"; backend=OnDisk())

    solenoid = Field{Face, Face, Center}(fts_u.grid)
    potential = Field{Center, Center, Center}(fts_u.grid)
    
    ψ = Observable(nov(solenoid[:, :, 1]))
    φ = Observable(nov(potential[:, :, 1]))
    b = @lift nov(fts_b[Time($t)][:, :, 1]) ./ (sp.L_ml * sp.f^2)
    
    on(t) do t
        u = interior(fts_u[Time(t)], :, :, 1)
        v = interior(fts_v[Time(t)], :, :, 1)
        
        ψ_new, φ_new = helmholtzdecomposition(u, v, sp.Lx / sp.Nx)
        
        set!(solenoid, ψ_new)
        fill_halo_regions!(solenoid)

        set!(potential, φ_new)
        fill_halo_regions!(potential)

        ψ[] = nov(solenoid[:, :, 1]) ./ (sp.L_ml^2 * sp.f^2)
        φ[] = nov(potential[:, :, 1]) ./ (sp.L_ml^2 * sp.f^2)
    end
    
    title = @lift begin
        t_day = @sprintf "%.0f" ($t / (24*3600))
        ft_str = @sprintf "%.1f" sp.f * $t
        γ_str =  @sprintf "%.2f" (1 / sp.T_mix / sp.f)
        L"\text{Potentials} \quad \gamma_\text{mix} = %$γ_str \quad t = %$t_day \, \text{day}\quad ft = %$ft_str"
    end

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

    ax_ψ = Axis(fig[2, 1]; ax_kw...)
    ax_φ = Axis(fig[2, 2]; ax_kw...)
    ax_b = Axis(fig[2, 3]; ax_kw...)

    hideydecorations!(ax_φ; ticks=false)
    hideydecorations!(ax_b; ticks=false)

    ht_ψ = begin
        xs = nov(xnodes(solenoid; with_halos=true)) ./ sp.L
        ys = nov(ynodes(solenoid; with_halos=true)) ./ sp.L
        data = ψ
        colormap = to_colormap(:curl)
        
        contourf!(ax_ψ, xs, ys, data; colormap, levels=range(-2, 2, 20), extendhigh=:auto, extendlow=:auto)
    end

    ht_φ = begin
        xs = nov(xnodes(potential; with_halos=true)) ./ sp.L
        ys = nov(ynodes(potential; with_halos=true)) ./ sp.L
        data = φ
        colormap = to_colormap(:balance)
        
        contourf!(ax_φ, xs, ys, data; colormap, levels=range(-0.1, 0.1, 10), extendhigh=:auto, extendlow=:auto)
    end

    ht_b = begin
        xs = nov(xnodes(fts_b; with_halos=true)) ./ sp.L
        ys = nov(ynodes(fts_b; with_halos=true)) ./ sp.L
        data = b
        heatmap!(ax_b, xs, ys, data; colormap=reverse(to_colormap(:deep)), colorrange=(-16, 16))
    end

    Colorbar(fig[3, 1], ht_ψ; label=L"\psi / L_\text{ml}^2 f^2", vertical=false, flipaxis=false)
    Colorbar(fig[3, 2], ht_φ; label=L"\phi / L_\text{ml}^2 f^2", vertical=false, flipaxis=false)
    Colorbar(fig[3, 3], ht_b; label=L"b / L_\text{ml} f^2", vertical=false, flipaxis=false)

    colgap!(fig.layout, 40)
    prettyrecord(n, fig, filename, frames; record_kw...)

    return fig
end

function field_video(run_id, frames, filename;
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
    fts_b = FieldTimeSeries(SURFACE, "b"; backend=OnDisk())

    u = @lift nov(fts_u[Time($t)][:, :, 1]) .* 100
    v = @lift nov(fts_v[Time($t)][:, :, 1]) .* 100
    b = @lift nov(fts_b[Time($t)][:, :, 1]) .* 100
    
    title = @lift begin
        t_day = @sprintf "%.0f" ($t / (24*3600))
        ft_str = @sprintf "%.1f" sp.f * $t
        γ_str = run_label(run_id)
        L"\text{Submesoscale} \quad \gamma_\text{mix} = %$γ_str \quad t = %$t_day \, \text{day}\quad ft = %$ft_str"
    end

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

    ax_u = Axis(fig[2, 1]; ax_kw...)
    ax_v = Axis(fig[2, 2]; ax_kw...)
    ax_b = Axis(fig[2, 3]; ax_kw...)

    hideydecorations!(ax_u; ticks=false)
    hideydecorations!(ax_v; ticks=false)

    ht_u = begin
        xs = nov(xnodes(fts_u; with_halos=true)) ./ sp.L
        ys = nov(ynodes(fts_u; with_halos=true)) ./ sp.L
        data = u
        heatmap!(ax_u, xs, ys, data; colormap=:balance, colorrange=(-50, 50))
    end

    ht_v = begin
        xs = nov(xnodes(fts_v; with_halos=true)) ./ sp.L
        ys = nov(ynodes(fts_v; with_halos=true)) ./ sp.L
        data = v
        heatmap!(ax_v, xs, ys, data; colormap=:balance, colorrange=(-50, 50))
    end

    ht_b = begin
        xs = nov(xnodes(fts_b; with_halos=true)) ./ sp.L
        ys = nov(ynodes(fts_b; with_halos=true)) ./ sp.L
        data = b
        heatmap!(ax_b, xs, ys, data; colormap=:thermal, colorrange=(0, 0.05))
    end

    Colorbar(fig[3, 1], ht_u; label=L"u / \text{cm}\,\text{s}^{-1}", vertical=false, flipaxis=false)
    Colorbar(fig[3, 2], ht_v; label=L"v / \text{cm}\,\text{s}^{-1}", vertical=false, flipaxis=false)
    Colorbar(fig[3, 3], ht_b; label=L"b / \text{cm}\,\text{s}^{-2}", vertical=false, flipaxis=false)

    colgap!(fig.layout, 40)
    prettyrecord(n, fig, filename, frames; record_kw...)

    return fig
end
