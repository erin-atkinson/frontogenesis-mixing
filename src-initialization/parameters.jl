# parameters.jl

default_inputs = (;
    # stop time and timescales for init, save
    stop_time = 100e4, save_time = 5e4,
    # base length and timescales
    f = 1e-4, L = 10e3, H = 1000,
    # grid size
    Nx = 1024, Ny = 1024, Nz = 128,
    # mixed layer
    βH_ml = 0.01,
    # baroclinic growth rate
    βσ = 0.01,
    # domain size
    βx = 10, βy = 10
)

@inline function create_simulation_parameters(input_parameters=(; ))
    ip = merge(default_inputs, input_parameters)

    # Stratification from deformation radius
    N² = (ip.L / ip.H * ip.f)^2
    H_ml = ip.βH_ml * ip.H
    # Thermocline BL growth rate estimate
    σ = ip.βσ * ip.f
    
    # Shear necessary for that 
    S = σ * sqrt(N²) / 0.3ip.f
    M² = S * ip.f
    U = S * ip.H

    # Domain size
    Lx = ip.βx * ip.L
    Ly = ip.βy * ip.L

    Lz = ip.H

    op = (; N², S, M², U, σ, Lx, Ly, Lz, H_ml)
    
    return merge(ip, op)
end

@inline function create_simulation_parameters(; input_parameters...)
    create_simulation_parameters(input_parameters)
end
