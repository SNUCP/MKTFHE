uniform_binary(N::Int64, T::Type) = 
    rand(ChaCha20Stream(), [T(0), T(1)], N)

uniform_ternary(N::Int64, T::Type) = 
    rand(ChaCha20Stream(), [-T(1), T(0), T(1)], N)

function block_binary(d::Int64, ℓ::Int64, T::Type)
    vec = Vector{T}(undef, d * ℓ)
    @inbounds @simd for i = 1 : d * ℓ
        vec[i] = T(0)
    end

    rng = ChaCha20Stream()

    @inbounds @simd for i = 0 : d - 1
        idx = rand(rng, 0 : ℓ)
        if idx != 0
            vec[i * ℓ + idx] = T(1)
        end
    end
    vec
end

gaussian(σ::Float64) =
    σ * randn(ChaCha20Stream(), Float64)

gaussian(N::Int64, σ::Float64) =
    σ * randn(ChaCha20Stream(), Float64, N)

uniform_random32(N::Int64) =
    rand(ChaCha20Stream(), UInt32, N)

uniform_random64(N::Int64) =
    rand(ChaCha20Stream(), UInt64, N)
