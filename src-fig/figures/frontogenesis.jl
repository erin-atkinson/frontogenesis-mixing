function frontogenesistimeseries(run_ids, colors=Makie.wong_colors())
    fig = Figure(; fontsize=18, size=(1400, 600))
    ax_Fb = Axis(fig[1, 1];
        xlabel = L"f t",
        ylabel = L"\langle |\nabla b|^2\rangle / f^4"
    )
    hidexdecorations!(ax_Fb; ticks=false, grid=false)
    
    ax_Fc = Axis(fig[2, 1];
        xlabel = L"f t",
        ylabel = L"\langle|\nabla c|^2\rangle / f^4"
    )

    ax_δb = Axis(fig[1, 2];
        xlabel = L"f t",
        ylabel = L"-\delta_b / f"
    )
    hidexdecorations!(ax_δb; ticks=false, grid=false)
    
    ax_δc = Axis(fig[2, 2];
        xlabel = L"f t",
        ylabel = L"-\delta_c / f"
    )

    ax_γσb = Axis(fig[1, 3];
        xlabel = L"f t",
        ylabel = L"-\sigma_b / f"
    )
    hidexdecorations!(ax_γσb; ticks=false, grid=false)
    
    ax_γσc = Axis(fig[2, 3];
        xlabel = L"f t",
        ylabel = L"-\sigma_c / f"
    )

    N = length(run_ids)
    
    lns = map(1:N) do n
        run_id = run_ids[n]
        color = colors[n]
        
        STRAINTENSOR = joinpath(scratchpath, run_id, "STRAINTENSOR.jld2")
        iterations, times = iterations_times(STRAINTENSOR)
        sp = simulation_parameters(STRAINTENSOR)
        
        Fb = timeseries_of(mean, STRAINTENSOR, "Fb", iterations)
        Fc = timeseries_of(mean, STRAINTENSOR, "Fc", iterations) 
        
        δb = meanproduct(STRAINTENSOR, iterations, "δ", "Fb") ./ Fb ./ sp.f
        δc = meanproduct(STRAINTENSOR, iterations, "δ", "Fc") ./ Fc ./ sp.f

        γσb = meanproduct(STRAINTENSOR, iterations, "γb", "σ", "Fb") ./ Fb ./ sp.f
        γσc = meanproduct(STRAINTENSOR, iterations, "γc", "σ", "Fc") ./ Fc ./ sp.f

        i_max = argmax(filt(Fb, 3))
        i_growth = argmax(diff(filt(Fb, 3)) ./ diff(times))

        t_max = times[i_max]
        F_max = Fb[i_max]
        i_growth = findfirst(Fb .> F_max / 2)
        i_decay = findlast(Fb .> F_max / 2)
        
        t_growth = times[i_growth]
        t_decay = times[i_decay]

        println("$n, $run_id: t_max = $(t_max * sp.f) / f, t_growth = $(t_growth * sp.f) / f, t_decay = $(t_decay * sp.f) / f")

        lines!(ax_Fb, times .* sp.f, Fb ./ sp.f^4; color)
        lines!(ax_Fc, times .* sp.f, Fc ./ sp.f^4; color)

        scatter!(ax_Fb, t_max * sp.f, Fb[i_max] / sp.f^4; color, markersize=8)
        scatter!(ax_Fb, t_growth * sp.f, Fb[i_growth] / sp.f^4; color, markersize=8)
        i_decay < length(times) && scatter!(ax_Fb, t_decay * sp.f, Fb[i_decay] / sp.f^4; color, markersize=8)
        scatter!(ax_Fb, t_max * sp.f, Fb[i_max] / sp.f^4; color=:white, markersize=5)
        scatter!(ax_Fb, t_growth * sp.f, Fb[i_growth] / sp.f^4; color=:white, markersize=5)
        i_decay < length(times) && scatter!(ax_Fb, t_decay * sp.f, Fb[i_decay] / sp.f^4; color=:white, markersize=5)
        
        lines!(ax_δb, times .* sp.f, -filt(δb, 1); color)
        lines!(ax_δc, times .* sp.f, -filt(δc, 1); color)
        
        lines!(ax_γσb, times .* sp.f, -filt(γσb, 1); color)
        lines!(ax_γσc, times .* sp.f, -filt(γσc, 1); color)
    end
    
    Legend(fig[1:2, 4], lns, map(run_label, run_ids), legend_title)
    
    return fig
end

function meanproduct(filename, iterations, args...)
    map(iterations) do iteration
        mean(mapreduce(f->get_field(filename, f, iteration), .*, args))
    end
end

function frontogenesisspectra(run_ids, ts, colors)
    fig = Figure(; fontsize=18, size=(1000, 600))
    ax_divergence = Axis(fig[1, 1];
        title = L"\text{Divergence}",
        xscale = log10, 
        xreversed = true, 
        xlabel = L"\lambda / \text{m}",
    )
    ax_strain = Axis(fig[1, 2];
        title = L"\text{Strain}",
        xscale = log10, 
        xreversed = true, 
        xlabel = L"\lambda / \text{m}",
    )
    ax_shear = Axis(fig[2, 1];
        title = L"\text{Shear dispersion}",
        xscale = log10, 
        xreversed = true, 
        xlabel = L"\lambda / \text{m}",
    )
    ax_flux = Axis(fig[2, 2];
        title = L"\text{Flux}",
        xscale = log10, 
        xreversed = true, 
        xlabel = L"\lambda / \text{m}",
    )
    ax_mixing = Axis(fig[1, 3];
        title = L"\text{Mixing}",
        xscale = log10, 
        xreversed = true, 
        xlabel = L"\lambda / \text{m}",
    )
    #hidexdecorations!(ax_divergence; ticks=false, grid=false)
    #hidexdecorations!(ax_strain; ticks=false, grid=false)
    
    lns = map(run_ids, ts, colors) do run_id, t, color
        FRONTOGENESIS = joinpath(scratchpath, run_id, "FRONTOGENESIS.jld2")
        iterations, times = iterations_times(FRONTOGENESIS)
        F = timeseries_of(a->mean(a.^2), FRONTOGENESIS, "bx", iterations)
        F += timeseries_of(a->mean(a.^2), FRONTOGENESIS, "by", iterations)
        
        sp = simulation_parameters(FRONTOGENESIS)
        i = argmin(abs.(t / sp.f .- times))
        println(i, ": ", times[i] * sp.f)
        iteration = iterations[i]

        δbx = get_field(FRONTOGENESIS, "δbx", iteration)
        δby = get_field(FRONTOGENESIS, "δby", iteration)

        γσbx = get_field(FRONTOGENESIS, "γσbx", iteration)
        γσby = get_field(FRONTOGENESIS, "γσby", iteration)
        
        mixing_bx = get_field(FRONTOGENESIS, "mixing_bx", iteration)
        mixing_by = get_field(FRONTOGENESIS, "mixing_by", iteration)

        flux_bx = get_field(FRONTOGENESIS, "flux_bx", iteration)
        flux_by = get_field(FRONTOGENESIS, "flux_by", iteration)
        
        diffusion_bx = get_field(FRONTOGENESIS, "diffusion_bx", iteration)
        diffusion_by = get_field(FRONTOGENESIS, "diffusion_by", iteration)
        
        bx = get_field(FRONTOGENESIS, "bx", iteration)
        by = get_field(FRONTOGENESIS, "by", iteration)
    
        λs, bxδbx  = modedecomposition(δbx, bx, sp.Lx / sp.Nx)
        λs, byδby  = modedecomposition(δby, by, sp.Lx / sp.Nx)
        λs, bxγσbx  = modedecomposition(γσbx, bx, sp.Lx / sp.Nx)
        λs, byγσby  = modedecomposition(γσby, by, sp.Lx / sp.Nx)
        
        λs, bx_mixing_bx = modedecomposition(mixing_bx, bx, sp.Lx / sp.Nx)
        λs, by_mixing_by = modedecomposition(mixing_by, by, sp.Lx / sp.Nx)
        
        λs, bx_flux_bx = modedecomposition(flux_bx, bx, sp.Lx / sp.Nx)
        λs, by_flux_by = modedecomposition(flux_by, by, sp.Lx / sp.Nx)

        λs, bx_diffusion_bx = modedecomposition(diffusion_bx, bx, sp.Lx / sp.Nx)
        λs, by_diffusion_by = modedecomposition(diffusion_by, by, sp.Lx / sp.Nx)
        
        λs, bxbx = modedecomposition(bx, bx, sp.Lx/sp.Nx)
        λs, byby = modedecomposition(by, by, sp.Lx/sp.Nx)
        #A = (bxbx .+ byby)
        A = 1e-8
        divergence = -(bxδbx .+ byδby) ./ A
        strain = -(bxγσbx .+ byγσby) ./ A

        shear = (bx_mixing_bx .+ by_mixing_by) ./ A
        flux = (bx_flux_bx .+ by_flux_by) ./ A
        mixing = (bx_diffusion_bx .+ by_diffusion_by) ./ A
        
        lines!(ax_shear, λs, filt(shear, 3) ./ sp.f; color)
        lines!(ax_flux, λs, filt(flux, 3) ./ sp.f; color)
        lines!(ax_mixing, λs, filt(mixing, 3) ./ sp.f; color)
        lines!(ax_divergence, λs, filt(divergence, 3) ./ sp.f; color)
        lines!(ax_strain, λs, filt(strain, 3) ./ sp.f; color)
    end
    
    Legend(fig[2, 3], lns, map(run_label, run_ids), legend_title; tellwidth=false)
    
    fig
end
