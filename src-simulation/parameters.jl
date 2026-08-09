# parameters.jl

default_inputs = (;
    # stop time and timescales for init, save
    stop_time = 100e4, spinup_time = 10e4, save_time = 5e4,
    # base length and timescales
    f = 1e-4, L = 10e3, H = 1000,
    # grid size
    Nx = 1024, Ny = 1024, Nz = 128,
    # mixed layer properties
    βL_ml = 0.3, βH_ml = 0.3,
    # baroclinic growth rate
    βσ = 0.1, 
    # surface forcing
    βB = 0.0, βτ = 0.0, θτ = 0.0,
    # domain size
    βx = 10, βy = 10
)

@inline function create_simulation_parameters(input_parameters=(; ))
    ip = merge(default_inputs, input_parameters)

    # Stratification from deformation radius
    N² = (ip.L / ip.H * ip.f)^2

    # Mixed layer
    L_ml = ip.βL_ml * ip.L
    H_ml = ip.βH_ml * ip.H
    N²_ml = (L_ml / H_ml * ip.f)^2
    
    H_th = ip.H - H_ml
    
    # Thermocline BL growth rate estimate
    σ = ip.βσ * ip.f
    
    # Shear necessary for that 
    S = σ * sqrt(N²) / 0.3ip.f
    M² = S * ip.f
    U = S * ip.H

    # Mixing rate
    B = ip.βB * H_ml^2 * ip.f^3
    T_mix = (H_ml^2 / B)^(1/3)
    ν = (H_ml^2 / T_mix) / 100
    
    τ = ip.βτ * ip.L^2 * ip.f^2
    
    # Thermodynamics...
    αV = 2.0678e-4 # K⁻¹
    cₚ = 4.1819e3 # J kg⁻¹ K⁻¹
    ρ = 1.027e3 # kg m⁻³
    g = 9.81 # m s⁻²
    Q = (cₚ * ρ) * B / (αV * g)

    # Domain size
    Lx = ip.βx * ip.L
    Ly = ip.βy * ip.L

    Lz = ip.H

    λ = 0.01

    start_time = -ip.spinup_time

    op = (; N², H_ml, L_ml, N²_ml, H_th, S, M², U, σ, B, ν, T_mix, τ, Q, Lx, Ly, Lz, λ, start_time)
    
    return merge(ip, op)
end

@inline function create_simulation_parameters(; input_parameters...)
    create_simulation_parameters(input_parameters)
end
