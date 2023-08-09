struct LWEkey{T<:Unsigned}
    n::Int64
    key::Vector{T}
    
    LWEkey(n::Int64, key::Vector{T}) where {T<:Unsigned} =
        new{T}(n, key)
end

LEVkey = LWEkey
GSWkey = LWEkey

binary_lwekey(n::Int64, T::Type) = 
    LWEkey(n, uniform_binary(n, T))

ternary_lwekey(n::Int64, T::Type) = 
    LWEkey(n, uniform_ternary(n, T))

block_binary_lwekey(d::Int64, ℓ::Int64, T::Type) =
    LWEkey(d * ℓ, block_binary(d, ℓ, T))

struct RLWEkey{T, R}
    k::Int64
    N::Int64
    key::Vector{<:RingPoly{T}}
    tkey::Vector{<:TransPoly{R}}
end

RLEVkey = RLWEkey
RGSWkey = RLWEkey
Unikey = RLWEkey

binary_ringkey(k::Int64, N::Int64, T::Type, ffter::FFTransformer{R}) where R = begin
    key = Vector{NativePoly{T}}(undef, k)
    tkey = Vector{TransNativePoly{R}}(undef, k)
    @inbounds @simd for i = 1 : k
        key[i] = NativePoly{T}(uniform_binary(N, T), N)
        tkey[i] = fft(key[i], ffter)
    end
    RLWEkey{T, R}(k, N, key, tkey)
end

ternary_ringkey(k::Int64, N::Int64, T::Type, ffter::FFTransformer{R}) where R = begin
    key = Vector{NativePoly{T}}(undef, k)
    tkey = Vector{TransNativePoly{R}}(undef, k)
    @inbounds @simd for i = 1 : k
        key[i] = NativePoly{T}(uniform_ternary(N, T), N)
        tkey[i] = fft(key[i], ffter)
    end
    RLWEkey{T, R}(k, N, key, tkey)
end

function partial_ringkey(k::Int64, N::Int64, lwekey::LWEkey{T}, ffter::FFTransformer{R}) where {T, R}
    n = lwekey.n
    key = Vector{NativePoly{T}}(undef, k)
    tkey = Vector{TransNativePoly{R}}(undef, k)
    @inbounds @simd for i = 1 : k
        if n ≥ N
            key[i] = NativePoly{T}(lwekey.key[(i-1)*N+1:i*N], N)
            n -= N
        elseif n ≥ 0
            key[i] = NativePoly{T}(vcat(lwekey.key[(i-1)*N+1:end], uniform_binary(N-n, T)), N)
            n -= N
        else
            key[i] = NativePoly{T}(uniform_binary(N, T), N)
        end
        tkey[i] = fft(key[i], ffter)
    end
    RLWEkey{T, R}(k, N, key, tkey)
end

function partial_ringkey(k::Int64, N::Int64, S::Type, lwekey::LWEkey{T}, ffter::FFTransformer{R}) where {T, R}
    n = lwekey.n
    key = Vector{NativePoly{S}}(undef, k)
    tkey = Vector{TransNativePoly{R}}(undef, k)
    @inbounds @simd for i = 1 : k
        if n ≥ N
            key[i] = NativePoly{S}(S.(lwekey.key[(i-1)*N+1:i*N]), N)
            n -= N
        elseif n ≥ 0
            key[i] = NativePoly{S}(vcat(S.(lwekey.key[(i-1)*N+1:end]), uniform_binary(N-n, S)), N)
            n -= N
        else
            key[i] = NativePoly{S}(uniform_binary(N, S), N)
        end
        tkey[i] = fft(key[i], ffter)
    end
    RLWEkey{S, R}(k, N, key, tkey)
end