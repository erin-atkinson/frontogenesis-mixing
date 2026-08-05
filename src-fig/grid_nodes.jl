include("../src-simulation/grid_faces.jl")

get_grid(file) = file["serialized/grid"]

function serialized_grid_nodes(file; with_halos=false, reshape=false)
    grid = get_grid(file)

    xsᶜ, ysᶜ, zsᶜ = nodes(grid, Center(), Center(), Center(); reshape, with_halos)
    xsᶠ, ysᶠ, zsᶠ = nodes(grid, Face(), Face(), Face(); reshape, with_halos)

    return xsᶜ, xsᶠ, ysᶜ, ysᶠ, zsᶜ, zsᶠ
end

nov = no_offset_view
@inline function halos(file)
    grid = file["serialized/grid"]
    return grid.Hx, grid.Hy, grid.Hz
end

function grid_nodes(file; with_halos=false, reshape=false)
    return serialized_grid_nodes(file; with_halos, reshape)
    
    xsᶜ = nov(file["grid/xᶜᵃᵃ"])
    xsᶠ = nov(file["grid/xᶠᵃᵃ"])

    ysᶜ = nov(file["grid/yᵃᶜᵃ"])
    ysᶠ = nov(file["grid/yᵃᶠᵃ"])

    zsᶠ = nov(file["grid/z/cᵃᵃᶠ"])
    zsᶜ = nov(file["grid/z/cᵃᵃᶜ"])

    if !with_halos
        Hx, Hy, Hz = halos(file)

        xsᶜ = xsᶜ[(Hx+1):(end-Hx)]
        xsᶠ = xsᶠ[(Hx+1):(end-Hx)]
        ysᶜ = ysᶜ[(Hy+1):(end-Hy)]
        ysᶠ = ysᶠ[(Hy+1):(end-Hy)]
        zsᶜ = zsᶜ[(Hz+1):(end-Hz)]
        zsᶠ = zsᶠ[(Hz+1):(end-Hz)]

        return xsᶜ, xsᶠ, ysᶜ, ysᶠ, zsᶜ, zsᶠ
    end

    xsᶜ, xsᶠ, ysᶜ, ysᶠ, zsᶜ, zsᶠ
end

grid_nodes(filename::String; kwargs...) = jldopen(file->grid_nodes(file; kwargs...), filename)
