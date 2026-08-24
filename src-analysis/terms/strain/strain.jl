@inline function shearstrain(i, j, k, grid, u, v)
    return ∂xᶠᶠᶜ(i, j, k, grid, v) + ∂yᶠᶠᶜ(i, j, k, grid, u)
end

@inline function normalstrain(i, j, k, grid, u, v)
    return ∂xᶜᶜᶜ(i, j, k, grid, u) - ∂yᶜᶜᶜ(i, j, k, grid, v)
end

# It's easiest to calculate this without doing the rotation explicitly
# I'm not going to worry too much about the grid error
@inline function strainefficiency(i, j, k, grid, u, v, b)
    bxᶜᶜᶜ = ℑxᶜᵃᵃ(i, j, k, grid, ∂xᶠᶜᶜ, b)
    bxᶠᶠᶜ = ℑyᵃᶠᵃ(i, j, k, grid, ∂xᶠᶜᶜ, b)
    byᶜᶜᶜ = ℑyᵃᶜᵃ(i, j, k, grid, ∂yᶜᶠᶜ, b)
    byᶠᶠᶜ = ℑxᶠᵃᵃ(i, j, k, grid, ∂yᶜᶠᶜ, b)

    bxbx = bxᶜᶜᶜ^2
    byby = byᶜᶜᶜ^2
    bxby = bxᶠᶠᶜ * byᶠᶠᶜ
    
    δ = ∂xᶜᶜᶜ(i, j, k, grid, u) + ∂yᶜᶜᶜ(i, j, k, grid, v)
    σₛ = ℑxyᶜᶜᵃ(i, j, k, grid, shearstrain, u, v)
    σₙ = normalstrain(i, j, k, grid, u, v)
    
    ux = ∂xᶜᶜᶜ(i, j, k, grid, u)
    vy = ∂yᶜᶜᶜ(i, j, k, grid, v)
    
    uy = ∂yᶠᶠᶜ(i, j, k, grid, u)
    vx = ∂xᶠᶠᶜ(i, j, k, grid, v)

    bxuxbx = ux * bxbx
    byuybx = ℑxyᶜᶜᵃ(i, j, k, grid, fg, uy, bxby)
    bxvxby = ℑxyᶜᶜᵃ(i, j, k, grid, fg, vx, bxby)
    byvyby = vy * byby
    
    biuijbj = bxuxbx + byuybx + bxvxby + byvyby
    γσ = biuijbj / (bxbx + byby) - δ
    σ = sqrt(σₛ^2 + σₙ^2)
    
    return γσ / σ
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

function StrainEfficiency(u, v, b)
    grid = u.grid
    loc = locationornothing((Center, Center, Center), u)
    return KernelFunctionOperation{loc...}(strainefficiency, grid, u, v, b)
end
