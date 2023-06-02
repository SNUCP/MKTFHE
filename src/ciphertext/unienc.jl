Uniparams{T} = GSWparams{T}
Uniparams_digit{T} = GSWparams_digit{T}

function decomptoith!(avec::Array{NativePoly{T}, 2}, a::Vector{NativePoly{T}}, i::Int64, params::Union{GSWparams_digit{T}, LEVparams_digit{T}}) where T
    @inbounds @simd for idx = 1 : i
        @inbounds @simd for j1 = 1 : a[1].N
            aj1 = divbits(a[idx].coeffs[j1], params.gveclog[end])
            @inbounds for j2 = params.l : -1 : 1
                tmp = aj1 & params.mask
                aj1 >>= params.logB
                carry = tmp >> (params.logB - 1)
                aj1 += carry
                tmp -= carry << params.logB
                avec[j2, idx].coeffs[j1] = tmp
            end
        end
    end
end

CRS{T} = Vector{<:RingPoly{T}}
TransCRS{T} = Vector{<:TransPoly{T}}

fft(a::CRS, ffter::FFTransformer{T}) where T = 
    fft.(a, Ref(ffter))

# UniEnc is not mutable, since you cannot generate a valid uni-encryption by homomorphic operations.
struct UniEnc{T<:Unsigned}
    d::Vector{RingPoly{T}}
    f::RLEV{T}

    function UniEnc(d::Vector{<:RingPoly{T}}, f::RLEV{T}) where {T<:Unsigned}
        new{T}(d, f)
    end
end

function unienc_encrypt(ta::TransCRS{S}, m::T, key::Unikey{T}, σ::Float64, ffter::FFTransformer{S}, params::Uniparams{T}) where {T, S}
    N = key.N
    r = ternary_ringkey(1, N, T, ffter)

    d = Vector{NativePoly{T}}(undef, params.l)
    tdi = bufftransnativepoly(N, S)

    @inbounds for i = 1 : params.l
        multo!(tdi, ta[i], r.tkey[1])
        d[i] = ifft(tdi, ffter)
        d[i].coeffs[1] += m * params.gvec[i]
        @inbounds @simd for j = 1 : key.N
            d[i].coeffs[j] += unsigned(signed(T)(round.(gaussian(σ))))
        end
    end

    f = rlev_encrypt(r.key[1], key, σ, ffter, params)

    UniEnc(d, f)
end

function unienc_encrypt(ta::TransCRS{S}, m::NativePoly{T}, key::Unikey{T}, σ::Float64, ffter::FFTransformer{S}, params::Uniparams{T}) where {T, S}
    N = key.N
    r = ternary_ringkey(1, N, T, ffter)

    d = Vector{NativePoly{T}}(undef, params.l)
    tdi = bufftransnativepoly(N, S)

    @inbounds for i = 1 : params.l
        multo!(tdi, ta[i], r.tkey[1])
        d[i] = ifft(tdi, ffter)
        @inbounds @simd for j = 1 : key.N
            d[i].coeffs[j] += m.coeffs[j] * params.gvec[i] + unsigned(signed(T)(round.(gaussian(σ))))
        end
    end

    f = rlev_encrypt(r.key[1], key, σ, ffter, params)

    UniEnc(d, f)
end

function gen_b(ta::TransCRS{S}, key::Unikey{T}, σ::Float64, ffter::FFTransformer{S}, params::Uniparams{T}) where {T, S}
    N = key.N
    b = Vector{NativePoly{T}}(undef, params.l)
    tbi = bufftransnativepoly(N, S)
    @inbounds @simd for i = 1 : params.l
        initialise!(tbi)
        mulsubto!(tbi, key.tkey[1], ta[i])
        b[i] = ifft(tbi, ffter)
        @inbounds @simd for j = 1 : N
            b[i].coeffs[j] += unsigned(signed(T)(round(gaussian(σ))))
        end
    end
    b
end

struct TransUniEnc{T}
    d::Vector{TransPoly{T}}
    f::TransRLEV{T}

    function TransUniEnc(d::Vector{<:TransPoly{T}}, f::TransRLEV{T}) where T
        new{T}(d, f)
    end
end

fft(ct::UniEnc, ffter::FFTransformer) = 
    TransUniEnc(fft.(ct.d, Ref(ffter)), fft(ct.f, ffter))

ifft(ct::TransUniEnc, ffter::FFTransformer) = 
    UniEnc(ifft.(ct.d, Ref(ffter)), ifft.(ct.f, ffter))