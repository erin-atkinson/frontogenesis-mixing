# Get mean state using an along-front average
include("terms/terms.jl")

fields = (:u, :v, :w, :b, :ke)

mean_fields = NamedTuple()
for ξ in fields
    ξ_bar = Symbol(ξ, :_bar)
    @eval begin
        $ξ_bar = Field(Average(input_fields.$ξ; dims=2))
        mean_fields = (; mean_fields..., $ξ_bar)
    end
end

skip_update = filter(a->a ∉ fields, keys(input_fields))

dependency_fields = mean_fields
output_fields = dependency_fields
