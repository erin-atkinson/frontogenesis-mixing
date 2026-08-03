# parameters.jl

default_inputs = (;
    # stop time and timescales for init, save
    stop_time = 100e4, spinup_time = 10e4, save_time = 5e4,
    # base length and timescales
    f = 1e-4, L = 10e3, H = 1e3,
    # grid size
    Nx = 1024, Ny = 1024, Nz = 128,
    # mixed layer ratios
    βL_ml = 0.1, βH_ml = 0.1,
    # surface forcing
    βB = 0.0, βτ = 0.0, θτ = 0.0,
    # domain size
    βx = 10, βy = 10,
    comment = ""
)

@inline function create_simulation_parameters(input_parameters=(; ))
    ip = merge(default_inputs, input_parameters)

    # Stratification from deformation radius
    N² = ip.L / ip.H * ip.f^2

    # Mixed layer deformation radius assumes Ri=1
    H_ml = βH_ml * ip.H
    L_ml = βL_ml * ip.L

    # This gives shear
    S = L_ml / H_ml * ip.f^2
    M² = S * ip.f
    U = S * H_ml

    # Thermocline BL growth rate estimate
    σ = 0.3 * S * ip.f / sqrt(N²)

    # Cooling and mixing rate
    B = ip.βB * L_ml^2 * ip.f^3
    T_mix = (ip.H^2 / B)^(1/3)

    # Wind from turbulence scale?
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

    λ = H_ml / 100

    start_time = -ip.spinup_time

    op = (; N², H_ml, L_ml, S, M², U, σ, B, T_mix, τ, Q, Lx, Ly, Lz, λ, start_time)
    
    return merge(ip, op)
end

@inline function create_simulation_parameters(; input_parameters...)
    create_simulation_parameters(input_parameters)
end
