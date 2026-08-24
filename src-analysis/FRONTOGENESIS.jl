include("terms/terms.jl")

include("terms/advection/advection.jl")
include("terms/advection/operators.jl")

include("terms/mixedlayer/mld.jl")
include("terms/strain/divergence.jl")
include("terms/strain/strain.jl")

mld = Field{Nothing, Nothing, Nothing}(grid)
set!(mld, sp.H_ml / 20)
fields = (:u, :v, :w, :b)

u = input_fields.u
v = input_fields.v
w = input_fields.w
b = input_fields.b

u_mld = Field(ML_Average(u, mld))
v_mld = Field(ML_Average(v, mld))
w_mld = Field(ML_Average(w, mld))
b_mld = Field(ML_Average(b, mld))

b_surface = Field(TopSlice(b))

mld_fields = (; u_mld, v_mld, w_mld, b_mld, b_surface)

loc = (Center(), Center(), Nothing())

mld_flux_density_x = Field(UcFlux(centered, u_mld, b_mld))
mld_flux_density_y = Field(VcFlux(centered, v_mld, b_mld))
mld_flux_density_z = Field(WcFlux(centered, w_mld, b_mld))
mld_flux_density = (; mld_flux_density_x, mld_flux_density_y, mld_flux_density_z)

advection_x = Field(@at loc -u_mld * ∂x(b_mld))
advection_y = Field(@at loc -v_mld * ∂y(b_mld))
advection = (; advection_x, advection_y)

flux_density_x = Field(ML_Average(UcFlux(centered, u, b), mld))
flux_density_y = Field(ML_Average(VcFlux(centered, v, b), mld))
flux_density_z = Field(ML_Average(WcFlux(centered, w, b), mld))
flux_density = (; flux_density_x, flux_density_y, flux_density_z)

turbulent_flux_density_x = Field(flux_density_x - mld_flux_density_x)
turbulent_flux_density_y = Field(flux_density_y - mld_flux_density_y)
turbulent_flux_density_z = Field(flux_density_z - mld_flux_density_z)
turbulent_flux_density = (; turbulent_flux_density_x, turbulent_flux_density_y, turbulent_flux_density_z)

mixing_x = Field(-∂x(turbulent_flux_density_x))
mixing_y = Field(-∂y(turbulent_flux_density_y))
mixing_z = Field(turbulent_flux_density_z / mld)
mixing = (; mixing_x, mixing_y, mixing_z)

diffusion = Field(-sp.ν * (b_surface - b_mld) / mld^2)
base = (; diffusion)

dependency_fields = merge(mld_fields, mld_flux_density, advection, flux_density, turbulent_flux_density, mixing, base)

# Now frontogenesis...
# I'm not going to worry about the position here

γ = Field(StrainEfficiency(u_mld, v_mld, b_mld))
δ = Field(HorizontalDivergence(u_mld, v_mld))
σₙ = Field(NormalStrain(u_mld, v_mld))
σₛ = Field(ShearStrain(u_mld, v_mld))
σ = Field(sqrt(σₙ^2 + σₛ^2))

bx = Field(∂x(b_mld))
by = Field(∂y(b_mld))

strain = (; γ, δ, σₙ, σₛ, σ, bx, by)

δbx = Field(bx * δ)
δby = Field(by * δ)
γσbx = Field(bx * γ * σ)
γσby = Field(by * γ * σ)

mixing_bx = Field(∂x(mixing_x) + ∂x(mixing_y))
mixing_by = Field(∂y(mixing_x) + ∂y(mixing_y))

flux_bx = Field(∂x(mixing_z))
flux_by = Field(∂y(mixing_z))

diffusion_bx = Field(∂x(diffusion))
diffusion_by = Field(∂y(diffusion))

frontogenesis = (; δbx, δby, γσbx, γσby, mixing_bx, mixing_by, flux_bx, flux_by, diffusion_bx, diffusion_by)
dependency_fields = merge(dependency_fields, strain, frontogenesis)

output_fields = (; frontogenesis..., bx, by, δ, γ, σ)

skip_update = filter(a->a ∉ fields, keys(input_fields))
