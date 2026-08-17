# Integrates over θ by binning
# c(k₁, k₂) dk₁ dk₂ = ξ(k) dk
function binspectrum(k₁, k₂, c)
    N_bins = length(k₁)÷2

    ks = range(0, maximum(k₁), N_bins + 1)

    ξs = zeros(ComplexF64, N_bins)

    ks2D = [sqrt(k^2 + l^2) for k in k₁, l in k₂]
    
    i_sorted = sortperm(ks2D[:])

    i = 1
    j = 1
    while i <= N_bins
        if ks2D[i_sorted[j]] > ks[i]
            i += 1
        else
            ξs[i] += c[i_sorted[j]]
            j += 1
        end
    end
    
    return (ks[1:end-1] + ks[2:end]) ./ 2, ξs
end

function modedecomposition(fA, A, Δx)
    c = conj.(fft(A)) .* fft(fA)
    k = fftfreq(length(A), 1/Δx)
    ks, ξs = binspectrum(k, k, c)
    return ks, real.(ξs)
end