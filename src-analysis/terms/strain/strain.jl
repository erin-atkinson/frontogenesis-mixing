@inline function shearstrain(i, j, k, grid, u, v)
    return ∂xᶠᶠᶜ(i, j, k, grid, v) + ∂yᶠᶠᶜ(i, j, k, grid, u)
end

@inline function normalstrain(i, j, k, grid, u, v)
    return ∂xᶜᶜᶜ(i, j, k, grid, u) - ∂yᶜᶜᶜ(i, j, k, grid, v)
end

@inline function strainefficiency(i, j, k, grid, b)
    bx² = ℑxᶜᵃᵃ(i, j, k, grid, c², ∂xᶠᶜᶜ, b)
    by² = ℑyᵃᶜᵃ(i, j, k, grid, c², ∂yᶜᶠᶜ, b)
    return (bx² - by²) / (bx² + by²)
end

function ShearStrain(u, v)
    grid = u.grid
    loc = locationornothing((Face, Face, Center), u)
    return KernelFunctionOperation{loc...}(shearstrain, grid, u, v)
end

function NormalStrain(u, v)
    grid = u.grid
    loc = locationornothing((Center, Center, Center), u)
    return KernelFunctionOperation{loc...}(normalstrain, grid, u, v)
end

function StrainEfficiency(b)
    grid = b.grid
    loc = locationornothing((Center, Center, Center), b)
    return KernelFunctionOperation{loc...}(strainefficiency, grid, b)
end
