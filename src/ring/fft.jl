function bit_reverse!(μ::Vector{T}) where T
    j = 0
    n = length(μ)
    for i = 1:n-1
        bit = n >> 1
        while j ≥ bit
            j -= bit
            bit >>= 1
        end
        j += bit
        if i < j
            μ[i+1], μ[j+1] = μ[j+1], μ[i+1]
        end
    end
end

# For reasons unknown, making FFTransformer mutable and setting every value it holds constant reduces a lot of allocations. 
mutable struct FFTransformer{T}
    const N::Integer
    const Ψ::Vector{Complex{T}}
    const Ψinv::Vector{Complex{T}}
    const roots::Vector{Complex{T}}
    const rootsinv::Vector{Complex{T}}
    const mask::Union{UInt32, UInt64}

    function FFTransformer{T}(N::Integer, bits::Integer) where {T<:AbstractFloat}
        @assert bits == 32 || bits == 64

        mask = bits == 32 ? 0xffffffff : 0xffffffffffffffff

        halfN = N >> 1
        idx = collect(0 : halfN-1)
        Ψ = Complex{T}.(exp.((-im * big(π) / halfN) .* idx))
        Ψinv = Complex{T}.(exp.((im * big(π) / halfN) .* idx))

        bit_reverse!(Ψ)
        bit_reverse!(Ψinv)

        # We scale rootsinv instead of dividing by N, for simplicity.
        roots = Complex{T}.(exp.((im * big(π) / N) .* idx))
        rootinv = Complex{T}.(exp.((-im * big(π) / N) .* idx) / halfN)

        new(N, Ψ, Ψinv, roots, rootinv, mask)
    end
end

# Twisting Z[X]/(Xᴺ+1) to Z(i)[X]/(X^(N/2)+1).
function fft(p::NativePoly, ffter::FFTransformer{T}) where {T<:AbstractFloat}
    N = p.N; halfN = N >> 1
    a = (signed.(p.coeffs[1:halfN]) - 1im * signed.(p.coeffs[halfN+1:end])) .* ffter.roots
    fft!(a, ffter.Ψ)

    TransNativePoly(a, N)
end

# Twisting Z[X]/(Xᴺ+1) to Z(i)[X]/(X^(N/2)+1).
function fftto!(t::TransNativePoly, p::NativePoly, ffter::FFTransformer{T}) where {T<:AbstractFloat}
    halfN = p.N >> 1
    @inbounds @simd for i = 1 : halfN
        t.coeffs[i] = (signed(p.coeffs[i]) - im * signed(p.coeffs[halfN + i])) * ffter.roots[i]
    end
    fft!(t.coeffs, ffter.Ψ)
end

# Twisting Z[X]/(Xᴺ+1) to Z(i)[X]/(X^(N/2)+1).
function ifft(p::TransNativePoly, ffter::FFTransformer{T}) where {T<:AbstractFloat}
    ifft!(p.coeffs, ffter.Ψinv)
    @. p.coeffs *= ffter.rootsinv

    NativePoly(native.(vcat(real(p.coeffs), -imag(p.coeffs)), [ffter.mask]), p.N)
end

# Twisting Z[X]/(Xᴺ+1) to Z(i)[X]/(X^(N/2)+1).
function ifftto!(p::NativePoly, t::TransNativePoly, ffter::FFTransformer{T}) where {T<:AbstractFloat}
    N = p.N; halfN = N >> 1
    ifft!(t.coeffs, ffter.Ψinv)
    @. t.coeffs *= ffter.rootsinv

    @. p.coeffs[1:halfN] = native(real(t.coeffs), ffter.mask)
    @. p.coeffs[halfN+1:end] = native(-imag(t.coeffs), ffter.mask)
end

# Cooley-Tukey over C[X]/(Xᴺ+1).
# Based on https://eprint.iacr.org/2016/504
# Note : I implemented several FFT algorithms, but this negacyclic version without bit-reversal was the most performant.
function fft!(a::Vector{T}, Ψ::Vector{T}) where {T<:Number}
    N = length(a)

    m = 1
    logkp1 = trailing_zeros(N)
    k = N >> 1
    while logkp1 > 0
        @inbounds @simd for i = 0 : m-1
            j1 = i << logkp1 + 1; j2 = j1 + k - 1
            @inbounds @simd for j = j1 : j2
                t, u = a[j], Ψ[m+i+1] * a[j+k]
                a[j], a[j+k] = t + u, t - u
            end
        end
        m <<= 1; logkp1 -= 1; k >>= 1
    end
end

# Gentleman-Sande over C[X]/(Xᴺ+1).
# Based on https://eprint.iacr.org/2016/504
function ifft!(a::Vector{T}, Ψinv::Vector{T}) where {T<:Number}
    N = length(a)

    m = N >> 1
    logkp1 = 1
    k = 1
    while m > 0
        @inbounds @simd for i = 0 : m - 1
            j1 = i << logkp1 + 1; j2 = j1 + k - 1
            @inbounds @simd for j = j1 : j2
                t, u = a[j], a[j+k]
                a[j], a[j+k] = t + u, Ψinv[m+i+1] * (t - u)
            end
        end
        m >>= 1; logkp1 += 1; k <<= 1
    end
end