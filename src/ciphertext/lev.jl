abstract type DECparams{T<:Unsigned} end
abstract type LEVparams{T} <: DECparams{T} end

struct LEVparams_digit{T} <: LEVparams{T}
    l::Int64
    logB::Int64
    halfB::T
    gveclog::Vector{Int64}
    gvec::Vector{T}
    mask::T

    function LEVparams_digit{T}(l::Int64, logB::Int64) where T
        B = T(1) << logB
        maxbits = bits(T)
        gveclog = maxbits .- collect(1:l) * logB
        gvec = T(1) .<< gveclog
        mask = B - 1
        new{T}(l, logB, B >> 1, gveclog, gvec, mask)
    end
end

mutable struct LEV{T<:Unsigned}
    const l::Int64
    stack::Vector{LWE{T}}

    function LEV(stack::Vector{LWE{T}}) where {T<:Unsigned}
        new{T}(length(stack), stack)
    end
end

function lev_encrypt(m::T, key::LEVkey{T}, σ::Float64, params::DECparams) where T
    stack = Vector{LWE{T}}(undef, params.l)
    @inbounds @simd for i = 1 : params.l
        stack[i] = lwe_encrypt(m * params.gvec[i], key, σ)
    end
    LEV(stack)
end

function lev_ith_encrypt(m::T, i::Int64, key::LEVkey{T}, σ::Float64, params::DECparams) where T
    stack = Vector{LWE{T}}(undef, params.l)
    @inbounds @simd for j = 1 : params.l
        stack[j] = lwe_ith_encrypt(m * params.gvec[j], i, key, σ)
    end
    LEV(stack)
end

add(x::LEV{T}, y::LEV{T}) where T = 
    LEV((@. add(x.stack, y.stack)))

add!(x::LEV{T}, y::LEV{T}) where T =
    addto!(x, x, y)

addto!(res::LEV{T}, x::LEV{T}, y::LEV{T}) where T = begin
    @inbounds @simd for i = 1 : res.l
        addto!(res.stack[i], x.stack[i], y.stack[i])
    end
end

sub(x::LEV{T}, y::LEV{T}) where T = 
    LEV((@. sub(x.stack, y.stack)))

sub!(x::LEV{T}, y::LEV{T}) where T =
    subto!(x, x, y)

subto!(res::LEV{T}, x::LEV{T}, y::LEV{T}) where T = begin
    @inbounds @simd for i = 1 : res.l
        subto!(res.stack[i], x.stack[i], y.stack[i])
    end
end

initialise!(x::LEV{T}) where T = begin
    @inbounds @simd for i = 1 : x.l
        initialise!(x.stack[i])
    end
end

RLEVparams_digit{T} = LEVparams_digit{T}

mutable struct RLEV{T<:Unsigned}
    const l::Int64
    stack::Vector{RLWE{T}}

    function RLEV(stack::Vector{RLWE{T}}) where {T<:Unsigned}
        new{T}(length(stack), stack)
    end
end

function rlev_encrypt(m::T, key::RLEVkey{T}, σ::Float64, ffter::FFTransformer{S}, params::DECparams) where {T<:Unsigned, S<:AbstractFloat}
    stack = Vector{RLWE{T}}(undef, params.l)
    @inbounds @simd for i = 1 : params.l
        stack[i] = rlwe_encrypt(params.gvec[i] * m, key, σ, ffter)
    end
    RLEV(stack)
end

function rlev_ith_encrypt(m::T, i::Int64, key::RLEVkey{T}, σ::Float64, ffter::FFTransformer{S}, params::DECparams) where {T<:Unsigned, S<:AbstractFloat}
    stack = Vector{RLWE{T}}(undef, params.l)
    @inbounds @simd for j = 1 : params.l
        stack[j] = rlwe_ith_encrypt(params.gvec[j] * m, i, key, σ, ffter)
    end
    RLEV(stack)
end

rlev_encrypt(m::NativePoly{T}, key::RLEVkey{T}, σ::Float64, ffter::FFTransformer{S}, params::DECparams) where {T<:Unsigned, S<:AbstractFloat} =
    RLEV([rlwe_encrypt(params.gvec[i] * m, key, σ, ffter) for i = 1 : params.l])

rlev_ith_encrypt(m::NativePoly{T}, i::Int64, key::RLEVkey{T}, σ::Float64, ffter::FFTransformer{S}, params::DECparams) where {T<:Unsigned, S<:AbstractFloat} =
    RLEV([rlwe_ith_encrypt(params.gvec[j] * m, i, key, σ, ffter) for j = 1 : params.l])

add(x::RLEV{T}, y::RLEV{T}) where T = 
    RLEV((@. add(x.stack, y.stack)))

add!(x::RLEV{T}, y::RLEV{T}) where T =
    addto!(x, x, y)

addto!(res::RLEV{T}, x::RLEV{T}, y::RLEV{T}) where T = begin
    @inbounds @simd for i = 1 : res.l
        addto!(res.stack[i], x.stack[i], y.stack[i])
    end
end

sub(x::RLEV{T}, y::RLEV{T}) where T = 
    RLEV((@. sub(x.stack, y.stack)))

sub!(x::RLEV{T}, y::RLEV{T}) where T =
    subto!(x, x, y)

subto!(res::RLEV{T}, x::RLEV{T}, y::RLEV{T}) where T = begin
    @inbounds @simd for i = 1 : res.l
        subto!(res.stack[i], x.stack[i], y.stack[i])
    end
end

initialise!(x::RLEV{T}) where T = begin
    @inbounds @simd for i = 1 : x.l
        initialise!(x.stack[i])
    end
end

mutable struct TransRLEV{T}
    const l::Int64
    stack::Vector{TransRLWE{T}}

    function TransRLEV(stack::Vector{TransRLWE{T}}) where T
        new{T}(length(stack), stack)
    end
end

add(x::TransRLEV{T}, y::TransRLEV{T}) where T = 
    TransRLEV((@. add(x.stack, y.stack)))

add!(x::TransRLEV{T}, y::TransRLEV{T}) where T =
    addto!(x, x, y)

addto!(res::TransRLEV{T}, x::TransRLEV{T}, y::TransRLEV{T}) where T = begin
    @inbounds @simd for i = 1 : res.l
        addto!(res.stack[i], x.stack[i], y.stack[i])
    end
end

sub(x::TransRLEV{T}, y::TransRLEV{T}) where T = 
    TransRLEV((@. sub(x.stack, y.stack)))

sub!(x::TransRLEV{T}, y::TransRLEV{T}) where T =
    subto!(x, x, y)

subto!(res::TransRLEV{T}, x::TransRLEV{T}, y::TransRLEV{T}) where T = begin
    @inbounds @simd for i = 1 : res.l
        subto!(res.stack[i], x.stack[i], y.stack[i])
    end
end

mul(x::TransPoly{T}, ct::TransRLEV{T}) where T = begin
    TransRLEV([mul(x, ct.stack[i]) for i = eachindex(ct.stack)])
end

mul!(x::TransPoly{T}, ct::TransRLEV{T}) where T = 
    multo!(ct, x, ct)

multo!(res::TransRLEV{T}, x::TransPoly{T}, ct::TransRLEV{T}) where T = begin
    @inbounds @simd for i = 1 : res.l
        multo!(res.stack[i], x, ct.stack[i])
    end
end

initialise!(x::TransRLEV{T}) where T = begin
    @inbounds @simd for i = 1 : x.l
        initialise!(x.stack[i])
    end
end

fft(ct::RLEV, ffter::FFTransformer) =
    TransRLEV((fft.(ct.stack, Ref(ffter))))

fftto!(res::TransRLEV, ct::RLEV, ffter::FFTransformer) = begin
    @inbounds @simd for i = 1 : res.l
        fftto!(res.stack[i], ct.stack[i], ffter)
    end
end

ifft(ct::TransRLEV, ffter::FFTransformer) = begin
    RLEV((ifft.(ct.stack, Ref(ffter))))
end

ifftto!(res::RLEV, ct::TransRLEV, ffter::FFTransformer) = begin
    @inbounds @simd for i = 1 : res.l
        ifftto!(res.stack[i], ct.stack[i], ffter)
    end
end