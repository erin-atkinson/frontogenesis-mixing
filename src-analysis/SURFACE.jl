include("terms/terms.jl")

include("terms/advection/advection.jl")
include("terms/advection/operators.jl")

include("terms/mixedlayer/mld.jl")


mld = Field{Nothing, Nothing, Nothing}(grid)
set!(mld, -sp.H_ml / 10)

u_mld = Field(ML_Average(u, mld))
v_mld = Field(ML_Average(v, mld))
w_at_mld = Field(ML_Interpolate(w, mld))
b_mld = Field(ML_Average(b, mld))
dbdz_at_mld = Field(ML_Interpolate(∂b∂z, mld))
b_at_mld = Field(ML_Interpolate(b, mld))

u_dag = Field(u - u_mld)
v_dag = Field(v - v_mld)
b_dag = Field(b - b_mld)

loc = (Center(), Center(), Nothing())

flux_density_x = Field(UcFlux(centered, u_mld, b_mld))
flux_density_y = Field(VcFlux(centered, v_mld, b_mld))
flux_density = (; flux_density_x, flux_density_y)

advection_x = Field(@at loc -u_mld * ∂x(b_mld))
advection_y = Field(@at loc -v_mld * ∂y(b_mld))
advection = (; advection_x, advection_y)

turbulent_flux_density_x = Field(ML_Average(UcFlux(centered, u_dag, b_dag)))
turbulent_flux_density_y = Field(ML_Average(VcFlux(centered, v_dag, b_dag)))
turbulent_flux_density = (; turbulent_flux_density_x, turbulent_flux_density_y)

mixing_x = Field(-∂x(turbulent_flux_density_x))
mixing_y = Field(-∂y(turbulent_flux_density_y))
mixing = (; mixing_x, mixing_y)


w_at_mld * (b_mld - b_at_mld) + 

mld_flux_density_x = Field(UcFlux(centered, u_mld, b_mld))
mld_advection_x = Field(@at (loc[1], loc[2], Nothing()) -u_mld * ∂x(b_mld))
mld_advection_background = Field(ML_Average(flux_density_background, constant_mld))

shear_dispersion_flux = Field(flux_density - mld_flux_density)
mld_turbulent_flux = Field(ML_Average(turbulent_flux_density_x, constant_mld))

shear_dispersion = Field(-∂x(shear_dispersion_flux))
mixing = Field(-∂x(mld_turbulent_flux))


