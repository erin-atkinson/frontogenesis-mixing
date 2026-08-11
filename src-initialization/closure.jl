νₕ = 100 * (sp.Lx / sp.Nx)^2 * sp.σ
νᵥ = 0.01 * (sp.Lz / sp.Nz)^2 * sp.σ
ν = νₕ
closure = (
    HorizontalScalarDiffusivity(; ν=νₕ, κ=νₕ),
    VerticalScalarDiffusivity(; ν=νᵥ, κ=νᵥ),
)
@info "Created closure"
