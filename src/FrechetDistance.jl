module FrechetDistance

using Distances
using LinearAlgebra

export frechet

"""
    frechet(P, Q)

Compute discrete Fréchet distance between polygonal curves `P` and `Q`.

`P` and `Q` are repesented as matrices with each column corresponding to a
point.

Adapts the algorithm from Eiter and Mannila, 1994,
*Computing Discrete Fréchet Distance*.
"""
function frechet(metric::PreMetric, P::AbstractMatrix, Q::AbstractMatrix)
    dP, m = size(P)
    dQ, n = size(Q)

    dP != dQ && throw(DimensionMismatch(
        "Points in polygonal lines P and Q must have the same dimension."))

    couplings = pairwise(metric, P, Q, dims=2)

    @inbounds for i in 2:m
        couplings[i, 1] = max(couplings[i-1, 1], couplings[i, 1])
    end

    @inbounds for j in 2:n
        couplings[1, j] = max(couplings[1, j-1], couplings[1, j])
        for i in 2:m
            carried_coupling = min(couplings[i-1, j-1],
                                   couplings[i, j-1],
                                   couplings[i-1, j])
            couplings[i, j] = max(carried_coupling, couplings[i, j])
        end
    end

    return couplings[m, n]
end

frechet(P, Q) = frechet(Euclidean(1e-12), P, Q)

end  # module