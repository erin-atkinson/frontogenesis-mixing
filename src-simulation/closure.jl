closure = if sp.ν ≈ 0
    nothing
else
    VerticalScalarDiffusivity(; sp.ν, κ=sp.ν)
end
@info "Created closure"
