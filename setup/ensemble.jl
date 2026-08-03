# Default
default_inputs = (;
    # stop time and timescales for init
    stop_time = 100e4, spinup_time = 10e4, save_time = 10e4,
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

profile_test_defaults = (;
    stop_time = 100e4, spinup_time = 10e4, save_time = 1e4,
    Nx = 64, Ny = 64,
    βx = 10 / 16, βy = 10 / 16
)

# Test profiles
profile_test_set = (;
    name = "PROFILE_TESTS",
    ips = [
        (; default_inputs..., profile_test_defaults..., βB = 0.01),
        (; default_inputs..., profile_test_defaults..., βB = 0.03),
        (; default_inputs..., profile_test_defaults..., βB = 0.1),
    ],
    filenames = [
        "profile-low-cooling",
        "profile-med-cooling",
        "profile-high-cooling",
    ]
)

ensemble = (; profile_test_set)
