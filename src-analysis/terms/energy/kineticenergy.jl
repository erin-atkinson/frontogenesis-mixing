@inline function kineticenergydensity(i, j, k, grid, u, v, w)
    u² = ℑxᶜᵃᵃ(i, j, k, grid, c², u)
    v² = ℑyᵃᶜᵃ(i, j, k, grid, c², v)
    w² = ℑzᵃᵃᶜ(i, j, k, grid, c², w)
    return (u² + v² + w²) / 2
end

function KineticEnergyDensity(u, v, w)
    grid = u.grid
    loc = locationornothing((Center, Center, Center), u)
    return KernelFunctionOperation{loc...}(kineticenergydensity, grid, u, v, w)
end
