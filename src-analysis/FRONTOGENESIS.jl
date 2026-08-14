include("terms/terms.jl")

include("terms/advection/advection.jl")
include("terms/advection/operators.jl")

include("terms/mixedlayer/mld.jl")

mld = Field{Nothing, Nothing, Nothing}(grid)
set!(mld, sp.H_ml / 10)
fields = (:u, :v, :w, :b)

u = input_fields.u
v = input_fields.v
w = input_fields.w
b = input_fields.b

u_mld = Field(ML_Average(u, mld))
v_mld = Field(ML_Average(v, mld))
b_mld = Field(ML_Average(b, mld))

∂b∂z = Field(∂z(b))
∂b∂z_at_mld = Field(ML_Interpolate(∂b∂z, mld))

b_at_mld = Field(ML_Interpolate(b, mld))
w_at_mld = Field(ML_Interpolate(w, mld))

mld_fields = (; u_mld, v_mld, b_mld, ∂b∂z, ∂b∂z_at_mld, b_at_mld, w_at_mld)

u_dag = Field(u - u_mld)
v_dag = Field(v - v_mld)
b_dag = Field(b - b_mld)
dag_fields = (; u_dag, v_dag, b_dag)

loc = (Center(), Center(), Nothing())

flux_density_x = Field(UcFlux(centered, u_mld, b_mld))
flux_density_y = Field(VcFlux(centered, v_mld, b_mld))
flux_density = (; flux_density_x, flux_density_y)

advection_x = Field(@at loc -u_mld * ∂x(b_mld))
advection_y = Field(@at loc -v_mld * ∂y(b_mld))
advection = (; advection_x, advection_y)

turbulent_flux_density_x = Field(ML_Average(UcFlux(centered, u_dag, b_dag), mld))
turbulent_flux_density_y = Field(ML_Average(VcFlux(centered, v_dag, b_dag), mld))
turbulent_flux_density = (; turbulent_flux_density_x, turbulent_flux_density_y)

mixing_x = Field(-∂x(turbulent_flux_density_x))
mixing_y = Field(-∂y(turbulent_flux_density_y))
mixing = (; mixing_x, mixing_y)

Fh = Field((w_at_mld * (b_mld - b_at_mld) + sp.ν * ∂b∂z_at_mld) / mld)
base = (; Fh)

dependency_fields = merge(mld_fields, dag_fields, flux_density, advection, turbulent_flux_density, mixing, base)

output_fields = (;
    u_mld,
    v_mld,
    b_mld,
    flux_density...,
    advection...,
    turbulent_flux_density...,
    mixing...,
    base...
)

skip_update = filter(a->a ∉ fields, keys(input_fields))
