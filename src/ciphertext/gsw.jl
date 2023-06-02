abstract type GSWparams{T} <: DECparams{T} end

struct GSWparams_digit{T} <: GSWparams{T}
    k::Int64
    l::Int64
    logB::Int64
    halfB::T
    gveclog::Vector{Int64}
    gvec::Vector{T}
    mask::T

    function GSWparams_digit{T}(k::Int64, l::Int64, logB::Int64) where T
        B = T(1) << logB
        maxbits = bits(T)
        gveclog = maxbits .- collect(1:l) * logB
        gvec = T(1) .<< gveclog
        mask = B - 1
        new{T}(k, l, logB, B >> 1, gveclog, gvec, mask)
    end
end

@inline decomp(a::T, params::LEVparams_digit{T}) where T = begin
    res = Vector{T}(undef, params.l)
    decompto!(res, a, params)
    res
end

@inline decomp(a::Vector{T}, params::GSWparams_digit{T}) where T = begin
    res = Array{T, 2}(undef, params.l, params.k)
    decompto!(res, a, params)
    res
end

@inline unbalanceddecompto!(avec::Vector{T}, a::T, params::Union{LEVparams_digit{T}, GSWparams_digit{T}}) where T = begin
    ai = divbits(a, params.gveclog[end])
    @inbounds for i = params.l : -1 : 1
        avec[i] = ai & params.mask
        ai >>= params.logB
    end
end

@inline decompto!(avec::Vector{T}, a::T, params::Union{LEVparams_digit{T}, GSWparams_digit{T}}) where T = begin
    ai = divbits(a, params.gveclog[end])
    @inbounds for i = params.l : -1 : 2
        avec[i] = ai & params.mask
        ai >>= params.logB
        ai += avec[i] >> (params.logB - 1)
        avec[i] -= (avec[i] & params.halfB) << 1
    end
    avec[1] = ai & params.mask
    avec[1] -= (avec[1] & params.halfB) << 1
end

@inline decompto!(avec::Array{T, 2}, a::Vector{T}, params::GSWparams_digit{T}) where T = begin
    ai = @. divbits(a, params.gveclog[end])
    @inbounds for i = params.l : -1 : 2
        @. avec[i, :] = ai & params.mask
        @. ai >>= params.logB
        @. ai += (avec[i, :] >> (params.logB - 1))
        @. avec[i, :] -= (avec[i, :] & params.halfB) << 1
    end
    @. avec[1, :] = ai & params.mask
    @. avec[1, :] -= (ai & params.halfB) << 1
end

@inline decomp(a::NativePoly{T}, params::Union{LEVparams_digit{T}, GSWparams_digit{T}}) where T = begin
    res = Vector{NativePoly{T}}(undef, params.l)
    @inbounds @simd for j = 1 : params.l
        res[j] = buffnativepoly(a.N, T)
    end
    decompto!(res, a, params)
    res
end

@inline decomp(a::Vector{NativePoly{T}}, params::GSWparams_digit{T}) where T = begin
    res = Array{NativePoly{T}}(undef, params.l, params.k)
    @inbounds @simd for i = 1 : params.k
        @inbounds @simd for j = 1 : params.l
            res[j, i] = buffnativepoly(a.N, T)
        end
    end
    decompto!(res, a, params)
    res
end

@inline decompto!(avec::Vector{NativePoly{T}}, a::NativePoly{T}, params::Union{LEVparams_digit{T}, GSWparams_digit{T}}) where T = begin
    @. avec[1].coeffs = divbits(a.coeffs, params.gveclog[end])
    @inbounds for j = params.l : -1 : 2
        @. avec[j].coeffs = avec[1].coeffs & params.mask
        @. avec[1].coeffs >>= params.logB
        @. avec[1].coeffs += avec[j].coeffs >> (params.logB - 1)
        @. avec[j].coeffs -= (avec[j].coeffs & params.halfB) << 1
    end
    @. avec[1].coeffs &= params.mask
    @. avec[1].coeffs -= (avec[1].coeffs & params.halfB) << 1
end

@inline decompto!(avec::Array{NativePoly{T}, 2}, a::Vector{NativePoly{T}}, params::GSWparams_digit{T}) where T = begin
    @inbounds @simd for idx = 1 : params.k
        @. avec[1, idx].coeffs = divbits(a[idx].coeffs, params.gveclog[end])
        @inbounds for j = params.l : -1 : 2
            @. avec[j, idx].coeffs = avec[1, idx].coeffs & params.mask
            @. avec[1, idx].coeffs >>= params.logB
            @. avec[1, idx].coeffs += avec[j, idx].coeffs >> (params.logB - 1)
            @. avec[j, idx].coeffs -= (avec[j, idx].coeffs & params.halfB) << 1
        end
        @. avec[1, idx].coeffs &= params.mask
        @. avec[1, idx].coeffs -= (avec[1, idx].coeffs & params.halfB) << 1
    end
end

mutable struct GSW{T<:Unsigned}
    const k::Int64
    basketb::LEV{T}
    basketa::Vector{LEV{T}}

    function GSW(basketb::LEV{T}, basketa::Vector{LEV{T}}) where {T<:Unsigned}
        new{T}(length(basketa), basketb, basketa)
    end
end

function gsw_encrypt(m::T, key::GSWkey{T}, σ::Float64, params::GSWparams) where T
    basketb = lev_encrypt(m, key, σ, params)
    basketa = Vector{LEV{T}}(undef, params.k)
    @inbounds @simd for i = 1 : params.k
        basketa[i] = lev_ith_encrypt(m, i, key, σ, params)
    end
    GSW(basketb, basketa)
end

add(x::GSW{T}, y::GSW{T}) where T = 
    GSW(add(x.basketb, y.basketb), (@. add(x.basket, y.basket)))

add!(x::GSW{T}, y::GSW{T}) where {T<:Unsigned} =
    addto!(x, x, y)

addto!(res::GSW{T}, x::GSW{T}, y::GSW{T}) where T = begin
    addto!(res.basketb, x.basketb, y.basketb)
    @inbounds @simd for i = eachindex(x.basket)
        addto!(res.basketa[i], x.basketa[i], y.basketa[i])
    end
end

sub(x::GSW{T}, y::GSW{T}) where T = 
    GSW(sub(x.basketb, y.basketb), (@. sub(x.basketa, y.basketa)))

sub!(x::GSW{T}, y::GSW{T}) where T =
    subto!(x, x, y)

subto!(res::GSW{T}, x::GSW{T}, y::GSW{T}) where T = begin
    subto!(res.basketb, x.basketb, y.basketb)
    @inbounds @simd for i = eachindex(x.basket)
        subto!(res.basketa[i], x.basketa[i], y.basketa[i])
    end
end

initialise!(x::GSW{T}) where T = begin
    initialise!(x.basketb)
    @inbounds @simd for i = 1 : x.k
        initialise!(x.basketa[i])
    end
end

mutable struct RGSW{T<:Unsigned}
    const k::Int64
    basketb::RLEV{T}
    basketa::Vector{RLEV{T}}

    function RGSW(basketb::RLEV{T}, basketa::Vector{RLEV{T}}) where {T<:Unsigned}
        new{T}(length(basketa), basketb, basketa)
    end
end

function rgsw_encrypt(m::T, key::RGSWkey{T}, σ::Float64, ffter::FFTransformer{S}, params::GSWparams) where {T, S}
    basketb = rlev_encrypt(m, key, σ, ffter, params)
    basketa = [rlev_ith_encrypt(m, i, key, σ, ffter, params) for i = 1 : params.k]
    RGSW(basketb, basketa)
end

function rgsw_encrypt(p::NativePoly{T}, key::RGSWkey{T}, σ::Float64, ffter::FFTransformer{S}, params::GSWparams) where {T, S}
    basketb = rlev_encrypt(p, key, σ, ffter, params)
    basketa = [rlev_ith_encrypt(p, i, key, σ, ffter, params) for i = 1 : params.k]
    RGSW(basketb, basketa)
end

add(x::RGSW{T}, y::RGSW{T}) where T = 
    RGSW(add(x.basketb, y.basketb), (@. add(x.basketa, y.basketa)))

add!(x::RGSW{T}, y::RGSW{T}) where T =
    addto!(x, x, y)

addto!(res::RGSW{T}, x::RGSW{T}, y::RGSW{T}) where T = begin
    addto!(res.basketb, x.basketb, y.basketb)
    @inbounds @simd for i = eachindex(x.basket)
        addto!(res.basketa[i], x.basketa[i], y.basketa[i])
    end
end

sub(x::RGSW{T}, y::RGSW{T}) where T = 
    RGSW(sub(x.basketb, y.basketb), (@. sub(x.basketa, y.basketa)))

sub!(x::RGSW{T}, y::RGSW{T}) where T =
    subto!(x, x, y)

subto!(res::RGSW{T}, x::RGSW{T}, y::RGSW{T}) where T = begin
    subto!(res.basketb, x.basketb, y.basketb)
    @inbounds @simd for i = eachindex(x.basket)
        subto!(res.basketa[i], x.basketa[i], y.basketa[i])
    end
end

initialise!(x::RGSW{T}) where T = begin
    initialise!(x.basketb)
    @inbounds @simd for i = 1 : x.k
        initialise!(x.basketa[i])
    end
end

mutable struct TransRGSW{T}
    const k::Int64
    basketb::TransRLEV{T}
    basketa::Vector{TransRLEV{T}}

    function TransRGSW(basketb::TransRLEV{T}, basketa::Vector{TransRLEV{T}}) where T
        new{T}(length(basketa), basketb, basketa)
    end 
end

add(x::TransRGSW{T}, y::TransRGSW{T}) where T = 
    TransRGSW(add(x.basketb, y.basketb), (@. add(x.basketa, y.basketa)))

add!(x::TransRGSW{T}, y::TransRGSW{T}) where T =
    addto!(x, x, y)

addto!(res::TransRGSW{T}, x::TransRGSW{T}, y::TransRGSW{T}) where T = begin
    addto!(res.basketb, x.basketb, y.basketb)
    @inbounds @simd for i = eachindex(x.basket)
        addto!(res.basketa[i], x.basketa[i], y.basketa[i])
    end
end

sub(x::TransRGSW{T}, y::TransRGSW{T}) where T = 
    TransRGSW(sub(x.basketb, y.basketb), (@. sub(x.basketa, y.basketa)))

sub!(x::TransRGSW{T}, y::TransRGSW{T}) where T =
    subto!(x, x, y)

subto!(res::TransRGSW{T}, x::TransRGSW{T}, y::TransRGSW{T}) where T = begin
    subto!(res.basketb, x.basketb, y.basketb)
    @inbounds @simd for i = eachindex(x.basket)
        subto!(res.basketa[i], x.basketa[i], y.basketa[i])
    end
end

initialise!(x::TransRGSW{T}) where T = begin
    initialise!(x.basketb)
    @inbounds @simd for i = 1 : x.k
        initialise!(x.basketa[i])
    end
end

fft(ct::RGSW, ffter::FFTransformer) =
    TransRGSW(fft(ct.basketb, ffter), (fft.(ct.basketa, Ref(ffter))))

fftto!(res::TransRGSW, ct::RGSW, ffter::FFTransformer) = begin
    fftto!(res.basketb, ct.basketb, ffter)
    @inbounds @simd for i = 1 : res.kp1
        fftto!(res.basketa[i], ct.basketa[i], ffter)
    end
end

ifft(ct::TransRGSW, ffter::FFTransformer) =
    RGSW(ifft(ct.basketb, ffter), (ifft.(ct.basketa, Ref(ffter))))

ifftto!(res::RGSW, ct::TransRGSW, ffter::FFTransformer) = begin
    ifftto!(res.basketb, ct.basketb, ffter)
    @inbounds @simd for i = 1 : res.kp1
        ifftto!(res.basketa[i], ct.basketa[i], ffter)
    end
end