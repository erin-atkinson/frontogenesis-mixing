# Default
default_inputs = (;
    # stop time and timescales for init
    stop_time = 0.0, spinup_time = 300e4, save_time = 10e4,
    # base length and timescales
    f = 1e-4, L = 50e3, H = 1000,
    # grid size
    Nx = 1024, Ny = 1024, Nz = 64,
    # mixed layer size
    βH_ml = 0.3, βL_ml = 0.1,
    # baroclinic growth rate
    βσ = 0.01,
    # surface forcing
    βB = 0.0, βτ = 0.0, θτ = 0.0,
    # domain size
    βx = 6, βy = 6,
    comment = ""
)

cooling_defaults = (;
    stop_time = 300e4, spinup_time = 0.0, save_time = 25e4,
)

# Test profiles
cooling_set = (;
    name = "MIXING",
    ips = [
        (; default_inputs..., cooling_defaults..., βB = 0.0, stop_time = 100e4),
        (; default_inputs..., cooling_defaults..., βB = 0.001, stop_time = 600e4),
        (; default_inputs..., cooling_defaults..., βB = 0.01, stop_time = 600e4),
        (; default_inputs..., cooling_defaults..., βB = 0.017, stop_time = 600e4),
        (; default_inputs..., cooling_defaults..., βB = 0.03, stop_time = 600e4),
        (; default_inputs..., cooling_defaults..., βB = 0.05, stop_time = 600e4),
        (; default_inputs..., cooling_defaults..., βB = 0.1, stop_time = 900e4),
    ],
    filenames = [
        "no-mixing",
        "tiny-mixing",
        "low-mixing",
        "medl-mixing",
        "med-mixing",
        "medh-mixing",
        "high-mixing",
    ]
)

ensemble = (; cooling_set)
init_dims = (; Nx=1024, Ny=1024, Nz=2, stop_time=1e7, spindown_time=0.0, save_time=2e5)
