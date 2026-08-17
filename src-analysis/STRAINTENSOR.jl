# Terms in the strain tensor for the surface velocity
include("terms/terms.jl")

include("terms/strain/vorticity.jl")
include("terms/strain/divergence.jl")
include("terms/strain/strain.jl")

fields = (:u, :v, :b)
frames = frames[1:5:end]

ζ = Field(VorticityZ(input_fields.u, input_fields.v))
δ = Field(HorizontalDivergence(input_fields.u, input_fields.v))
σₙ = Field(NormalStrain(input_fields.u, input_fields.v))
σₛ = Field(ShearStrain(input_fields.u, input_fields.v))
σ = Field(sqrt(σₙ^2 + σₛ^2))
γ = Field(StrainEfficiency(input_fields.b))

dependency_fields = (; ζ, δ, σₙ, σₛ, σ, γ)
output_fields = dependency_fields
skip_update = filter(a->a ∉ fields, keys(input_fields))
