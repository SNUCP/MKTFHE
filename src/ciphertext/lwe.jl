mutable struct LWE{T<:Unsigned}
    const n::Int64
    b::T
    a::Vector{T}

    function LWE(b::T, a::Vector{T}) where {T<:Unsigned}
        new{T}(length(a), b, a)
    end
end

LWEsample(key::LWEkey{T}, σ::Float64) where T = begin
    e = unsigned(round(signed(T), gaussian(σ)))
    a = rand(ChaCha20Stream(), T, key.n)
    b = T(-reduce(+, a .* key.key) + e)
    LWE(b, a)
end

lwe_encrypt(m::T, key::LWEkey{T}, σ::Float64) where T = begin
    res = LWEsample(key, σ)
    res.b += m
    res
end

lwe_ith_encrypt(m::T, i::Int64, key::LWEkey{T}, σ::Float64) where T = begin
    res = LWEsample(key, σ)
    res.a[i] += m
    res
end

phase(ct::LWE{T}, key::LWEkey{T}) where T = begin
    reduce(+, ct.a .* key.key) + ct.b
end

add(x::LWE{T}, y::LWE{T}) where T = 
    LWE(x.b + y.b, (@. x.a + y.a))

add!(x::LWE{T}, y::LWE{T}) where T =
    addto!(x, x, y)

addto!(res::LWE{T}, x::LWE{T}, y::LWE{T}) where T = begin
    res.b = x.b + y.b
    @. res.a = x.a + y.a
end

sub(x::LWE{T}, y::LWE{T}) where T = 
    LWE(x.b - y.b, (@. x.a - y.a))

sub!(x::LWE{T}, y::LWE{T}) where T = 
    subto!(x, x, y)

subto!(res::LWE{T}, x::LWE{T}, y::LWE{T}) where T = begin
    res.b = x.b - y.b
    @. res.a = x.a - y.a
end

initialise!(x::LWE{T}) where T = begin
    x.b = 0 * x.b
    @. x.a = 0 * x.a
end

mutable struct RLWE{T<:Unsigned}
    const k::Int64
    const N::Int64
    b::RingPoly{T}
    a::Vector{<:RingPoly{T}}

    function RLWE(b::S, a::Vector{S}) where {T<:Unsigned, S<:RingPoly{T}}
        @assert b.N == a[1].N
        new{T}(length(a), b.N, b, a)
    end

    function RLWE(b::S, a::S) where {T<:Unsigned, S<:RingPoly{T}}
        @assert b.N == a.N
        new{T}(1, b.N, b, [a])
    end
end

function RLWEsample(key::RLWEkey{T}, σ::Float64, ffter::FFTransformer{S}) where {T, S}
    N = key.N
    a = [randnativepoly(N, T) for _ = 1 : key.k]
    b = buffnativepoly(N, T)
    tb = zerotransnativepoly(N, S)
    ta = zerotransnativepoly(N, S)
    @inbounds for i = 1 : key.k
        fftto!(ta, a[i], ffter)
        mulsubto!(tb, key.tkey[i], ta)
    end
    ifftto!(b, tb, ffter)
    @inbounds @simd for i = 1 : N
        b.coeffs[i] += unsigned(round(signed(T), gaussian(σ)))
    end
    RLWE(b, a)
end

rlwe_encrypt(m::T, key::RLWEkey{T}, σ::Float64, ffter::FFTransformer{S}) where {T, S} = begin
    res = RLWEsample(key, σ, ffter)
    res.b.coeffs[1] += m
    res
end

rlwe_ith_encrypt(m::T, i::Int64, key::RLWEkey{T}, σ::Float64, ffter::FFTransformer{S}) where {T, S} = begin
    res = RLWEsample(key, σ, ffter)
    res.a[i].coeffs[1] += m
    res
end

rlwe_encrypt(p::NativePoly{T}, key::RLWEkey{T}, σ::Float64, ffter::FFTransformer{S}) where {T, S} = begin
    res = RLWEsample(key, σ, ffter)
    add!(res.b, p)
    res
end

rlwe_ith_encrypt(p::NativePoly{T}, i::Int64, key::RLWEkey{T}, σ::Float64, ffter::FFTransformer{S}) where {T, S} = begin
    res = RLWEsample(key, σ, ffter)
    add!(res.a[i], p)
    res
end

function phase(ct::RLWE{T}, key::RLWEkey{T}, ffter::FFTransformer{S}) where {T, S}
    N = ct.N
    tres = zerotransnativepoly(N, S)
    ta = zerotransnativepoly(N, S)
    tkey = [TransNativePoly(Complex{S}.(key.tkey[i].coeffs), N) for i = 1 : key.k]
    fftto!(tres, ct.b, ffter)
    @inbounds for i = 1 : key.k
        fftto!(ta, ct.a[i], ffter)
        muladdto!(tres, tkey[i], ta)
    end
    ifft(tres, ffter)
end

add(x::RLWE{T}, y::RLWE{T}) where T = 
    RLWE(add(x.b, y.b), (@. add(x.a, y.a)))

add!(x::RLWE{T}, y::RLWE{T}) where T =
    addto!(x, x, y)

addto!(res::RLWE{T}, x::RLWE{T}, y::RLWE{T}) where T = begin
    addto!(res.b, x.b, y.b)
    @inbounds @simd for i = 1 : res.k
        addto!(res.a[i], x.a[i], y.a[i])
    end
end

sub(x::RLWE{T}, y::RLWE{T}) where T = 
    RLWE(sub(x.b, y.b), (@. sub(x.a, y.a)))

sub!(x::RLWE{T}, y::RLWE{T}) where T =
    subto!(x, x, y)

subto!(res::RLWE{T}, x::RLWE{T}, y::RLWE{T}) where T = begin
    subto!(res.b, x.b, y.b)
    @inbounds @simd for i = 1 : res.k
        subto!(res.a[i], x.a[i], y.a[i])
    end
end

initialise!(x::RLWE{T}) where T = begin
    initialise!(x.b)
    @inbounds @simd for i = 1 : x.k
        initialise!(x.a[i])
    end
end

mutable struct TransRLWE{T}
    const k::Int64
    b::TransPoly{T}
    a::Vector{TransPoly{T}}

    function TransRLWE(b::S, a::Vector{S}) where {T, S<:TransPoly{T}}
        @assert a[1].N == b.N
        new{T}(length(a), b, a)
    end

    function TransRLWE(b::S, a::S) where {T, S<:TransPoly{T}}
        @assert a.N == b.N
        new{T}(1, b, [a])
    end
end

add(x::TransRLWE{T}, y::TransRLWE{T}) where T = 
    TransRLWE(add(x.b, y.b), (@. add(x.a, y.a)))

add!(x::TransRLWE{T}, y::TransRLWE{T}) where T =
    addto!(x, x, y)

addto!(res::TransRLWE{T}, x::TransRLWE{T}, y::TransRLWE{T}) where T = begin
    addto!(res.b, x.b, y.b)
    @inbounds @simd for i = 1 : res.k
        addto!(res.a[i], x.a[i], y.a[i])
    end
end

sub(x::TransRLWE{T}, y::TransRLWE{T}) where T = 
    TransRLWE(sub(x.b, y.b), (@. sub(x.a, y.a)))

sub!(x::TransRLWE{T}, y::TransRLWE{T}) where T =
    subto!(x, x, y)

subto!(res::TransRLWE{T}, x::TransRLWE{T}, y::TransRLWE{T}) where T = begin
    subto!(res.b, x.b, y.b)
    @inbounds @simd for i = 1 : res.k
        subto!(res.a[i], x.a[i], y.a[i])
    end
end

mul(x::TransPoly{T}, ct::TransRLWE{T}) where T = begin
    TransRLWE(mul(x, ct.b), [mul(x, ct.a[i]) for i = 1 : ct.k])
end

mul!(x::TransPoly{T}, ct::TransRLWE{T}) where T = 
    multo!(ct, x, ct)

multo!(res::TransRLWE{T}, x::TransPoly{T}, ct::TransRLWE{T}) where T = begin
    multo!(res.b, x, ct.b)
    @inbounds @simd for i = 1 : res.k
        multo!(res.a[i], x, ct.a[i])
    end
end

muladdto!(res::TransRLWE{T}, x::TransPoly{T}, ct::TransRLWE{T}) where T = begin
    muladdto!(res.b, x, ct.b)
    @inbounds @simd for i = 1 : res.k
        muladdto!(res.a[i], x, ct.a[i])
    end
end

mulsubto!(res::TransRLWE{T}, x::TransPoly{T}, ct::TransRLWE{T}) where T = begin
    mulsubto!(res.b, x, ct.b)
    @inbounds @simd for i = 1 : res.k
        mulsubto!(res.a[i], x, ct.a[i])
    end
end

initialise!(x::TransRLWE{T}) where T = begin
    initialise!(x.b)
    @inbounds @simd for i = 1 : x.k
        initialise!(x.a[i])
    end
end

fft(ct::RLWE, ffter::FFTransformer) =
    TransRLWE(fft(ct.b, ffter), (fft.(ct.a, Ref(ffter))))

fftto!(res::TransRLWE, ct::RLWE, ffter::FFTransformer) = begin
    fftto!(res.b, ct.b, ffter)
    @inbounds @simd for i = 1 : res.k
        fftto!(res.a[i], ct.a[i], ffter)
    end
end

ifft(ct::TransRLWE, ffter::FFTransformer) = 
    RLWE(ifft(ct.b, ffter), (ifft.(ct.a, Ref(ffter))))

ifftto!(res::RLWE, ct::TransRLWE, ffter::FFTransformer) = begin
    ifftto!(res.b, ct.b, ffter)
    @inbounds @simd for i = 1 : res.k
        ifftto!(res.a[i], ct.a[i], ffter)
    end
end