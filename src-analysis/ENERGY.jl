include("terms/terms.jl")
include("terms/energy/kineticenergy.jl")
include("terms/energy/potentialenergy.jl")

fields = (:u, :v, :w, :b)

ke_density = Field(KineticEnergyDensity(input_fields.u, input_fields.v, input_fields.w))
pe_density = Field(PotentialEnergyDensity(input_fields.b))

ke = Field(Integral(ke_density))
pe = Field(Integral(pe_density))

dependency_fields = (; ke_density, pe_density, ke, pe)
output_fields = dependency_fields
skip_update = filter(a->a ∉ fields, keys(input_fields))
