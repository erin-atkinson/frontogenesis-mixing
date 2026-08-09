# Default
default_inputs = (;
    # stop time and timescales for init
    stop_time = 0.0, spinup_time = 300e4, save_time = 10e4,
    # base length and timescales
    f = 1e-4, L = 50e3, H = 1000,
    # grid size
    Nx = 1024, Ny = 1024, Nz = 128,
    # mixed layer size
    βH_ml = 0.3, βL_ml = 0.1,
    # baroclinic growth rate
    βσ = 0.01,
    # surface forcing
    βB = 0.0, βτ = 0.0, θτ = 0.0,
    # domain size
    βx = 10, βy = 10,
    comment = ""
)

profile_test_defaults = (;
    stop_time = 100e4, spinup_time = 0.0, save_time = 1e4,
    Nx = 64, Ny = 64,
    βx = 10 / 16, βy = 10 / 16
)

# Test profiles
profile_test_set = (;
    name = "PROFILE_TESTS",
    ips = [
        (; default_inputs..., profile_test_defaults..., βB = 0.3),
        (; default_inputs..., profile_test_defaults..., βB = 1.0),
        (; default_inputs..., profile_test_defaults..., βB = 3.0),
    ],
    filenames = [
        "profile-low-cooling",
        "profile-med-cooling",
        "profile-high-cooling",
    ]
)

cooling_defaults = (;
    stop_time = 100e4, spinup_time = 0.0, save_time = 10e4,
    Nx = 1024, Ny = 1024,
    βx = 10, βy = 10
)

# Test profiles
cooling_set = (;
    name = "COOLING",
    ips = [
        (; default_inputs..., cooling_defaults..., βB = 0.0),
        (; default_inputs..., cooling_defaults..., βB = 0.3),
        (; default_inputs..., cooling_defaults..., βB = 1.0),
        (; default_inputs..., cooling_defaults..., βB = 3.0),
    ],
    filenames = [
        "no-mixing",
        "low-mixing",
        "med-mixing",
        "high-mixing",
    ]
)

ensemble = (; profile_test_set, cooling_set)
init_dims = (; Nx=64, Ny=64, Nz=16, stop_time=6e7, spindown_time=1e6, save_time=1e6)
