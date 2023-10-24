# T is unsigned integer for (R)LWE, R, S are Float type for FFT operations in KeyGen and BlindRotate, respectively.
abstract type TFHEparams{T<:Unsigned, R<:AbstractFloat, S<:AbstractFloat} end

abstract type CGGIparams{T, R, S} <: TFHEparams{T, R, S} end

struct TFHEparams_bin{T, R, S} <: CGGIparams{T, R, S}
    n::Int64            # LWE dimension
    α::Float64          # LWE noise standard deviation
    
    f::Int64            # key-switching gadget length
    logD::Int64         # key-switching gadget size

    N::Int64            # RLWE dimension
    k::Int64            # RLWE length
    β::Float64          # RLWE noise standard deviation
    
    l_gsw::Int64        # blind-rotation gadget length
    logB_gsw::Int64     # blind-rotation gadget size
end

# LMSS23 : Faster TFHE Bootstrapping from Block Binary Distribution
struct TFHEparams_block{T, R, S} <: CGGIparams{T, R, S}
    d::Int64        # number of the blocks
    ℓ::Int64        # length of each block
    α::Float64      # LWE noise standard deviation
    
    f::Int64        # key-switching gadget length
    logD::Int64     # key-switching gadget size

    N::Int64        # RLWE dimension
    k::Int64        # RLWE length
    β::Float64      # RLWE noise standard deviation
    
    l_gsw::Int64    # blind-rotation gadget length
    logB_gsw::Int64 # blind-rotation gadget size
end

abstract type MKTFHEparams{T, R, S} <: TFHEparams{T, R, S} end

struct CCSparams{T, R, S} <: MKTFHEparams{T, R, S}
    n::Int64        # LWE dimension
    α::Float64      # LWE noise standard deviation
    
    f::Int64        # key-switching gadget length
    logD::Int64     # key-switching gadget size

    N::Int64        # RLWE dimension
    β::Float64      # RLWE noise standard deviation
    
    l_uni::Int64    # blind-rotation gadget length
    logB_uni::Int64 # blind-rotation gadget size

    k::Int64        # number of the parties 
end

# T is unsigned for LWE, R is unsigned for RLWE and U, S are Float for FFT opertaions in KeyGen and BlindRotate, respectively.
struct KMSparams{T, R, U, S} <: MKTFHEparams{R, U, S} where {T<:Unsigned}
    n::Int64            # LWE dimension
    α::Float64          # LWE noise standard deviation
    
    f::Int64            # key-switching gadget length
    logD::Int64         # key-switching gadget size

    N::Int64            # RLWE dimension
    β::Float64          # RLWE noise standard deviation
    
    l_gsw::Int64        # GSW gadget length
    logB_gsw::Int64     # GSW gadget size

    l_lev::Int64        # LEV gadget length
    logB_lev::Int64     # LEV gadget size

    l_uni::Int64        # UniEnc gadget length
    logB_uni::Int64     # UniEnc gadget size

    k::Int64            # number of the parties 
end

# T is unsigned for LWE, R is unsigned for RLWE and U, S are Float for FFT opertaions in KeyGen and BlindRotate, respectively.
struct KMSparams_block{T, R, U, S} <: MKTFHEparams{R, U, S} where {T<:Unsigned}
    d::Int64        # number of the blocks
    ℓ::Int64        # length of each block
    α::Float64      # LWE noise standard deviation
    
    f::Int64        # key-switching gadget length
    logD::Int64     # key-switching gadget size

    N::Int64        # RLWE dimension
    β::Float64      # RLWE noise standard deviation
    
    l_gsw::Int64    # blind-rotation gadget length
    logB_gsw::Int64 # blind-rotation gadget size
    
    l_lev::Int64        # LEV gadget length
    logB_lev::Int64     # LEV gadget size

    l_uni::Int64        # UniEnc gadget length
    logB_uni::Int64     # UniEnc gadget size

    k::Int64            # number of the parties 
end

abstract type TFHEscheme{T<:Unsigned, S<:AbstractFloat} end

abstract type SKscheme{T, S} <: TFHEscheme{T, S} end

struct CGGI{T, S} <: SKscheme{T, S}
    k::Int64
    n::Int64
    N::Int64
    kskpar::LEVparams_digit{T}
    gswpar::GSWparams_digit{T}
    ffter::FFTransformer{S}
    monomial::Vector{TransNativePoly{S}}
    btk::BootKey_bin{T}
end

keygen_params(params::TFHEparams_bin{T, R, S}, ffter::FFTransformer{R}) where {T, R, S} =
    binary_lwekey(params.n, T), binary_ringkey(params.k, params.N, T, ffter)
    
function getmonomial(T::Type, ffter::FFTransformer{S}) where S
    N = ffter.N

    monomial = Vector{TransNativePoly{S}}(undef, 2N)
    monomial[2N] = zerotransnativepoly(N, S)

    tmppoly = zeronativepoly(N, T)
    tmppoly.coeffs[1] = -T(1)

    @inbounds for i = 2 : N
        tmppoly.coeffs[i] = T(1)
        monomial[i-1] = fft(tmppoly, ffter)
        tmppoly.coeffs[i] = T(0)
    end

    tmppoly.coeffs[1] = -T(2)
    monomial[N] = fft(tmppoly, ffter)
    tmppoly.coeffs[1] = -T(1)
    @inbounds for i = 2 : N
        tmppoly.coeffs[i] = -T(1)
        monomial[N+i-1] = fft(tmppoly, ffter)
        tmppoly.coeffs[i] = T(0)
    end

    monomial
end

"""
setup outputs LWE key, RLWE key, and Scheme.
"""
function setup(params::TFHEparams_bin{T, R, S}) where {T, R, S}
    k, n, N = params.k, params.n, params.N

    ffter = FFTransformer{S}(N, bits(T))
    ffter_keygen = FFTransformer{R}(N, bits(T))

    lwekey, ringkey = keygen_params(params, ffter_keygen)

    kskpar = LEVparams_digit{T}(params.f, params.logD)
    gswpar = GSWparams_digit{T}(k, params.l_gsw, params.logB_gsw)
    
    monomial = getmonomial(T, ffter)
    btk = BootKey_bin(lwekey, ringkey, kskpar, params.α, gswpar, params.β, ffter_keygen, ffter)
    
    lwekey, ringkey, CGGI{T, S}(k, n, N, kskpar, gswpar, ffter, monomial, btk)
end

struct LMSS{T, S} <: SKscheme{T, S}
    k::Int64
    ℓ::Int64
    d::Int64
    n::Int64
    N::Int64
    kskpar::LEVparams_digit{T}
    gswpar::GSWparams_digit{T}
    ffter::FFTransformer{S}
    monomial::Vector{TransNativePoly{S}}
    btk::BootKey_block{T}
end

keygen_params(params::TFHEparams_block{T, R, S}, ffter::FFTransformer{R}) where {T, R, S} = begin
    lwekey = block_binary_lwekey(params.d, params.ℓ, T)
    ringkey = partial_ringkey(params.k, params.N, lwekey, ffter)
    lwekey, ringkey
end

"""
setup outputs LWE key, RLWE key, and Scheme.
"""
function setup(params::TFHEparams_block{T, R, S}) where {T, R, S}
    k, n, N = params.k, params.ℓ * params.d, params.N

    ffter = FFTransformer{S}(N, bits(T))
    ffter_keygen = FFTransformer{R}(N, bits(T))

    lwekey, ringkey = keygen_params(params, ffter_keygen)

    kskpar = LEVparams_digit{T}(params.f, params.logD)
    gswpar = GSWparams_digit{T}(k, params.l_gsw, params.logB_gsw)
    
    monomial = getmonomial(T, ffter)
    btk = BootKey_block(lwekey, ringkey, kskpar, params.α, gswpar, params.β, ffter_keygen, ffter)
    
    lwekey, ringkey, LMSS{T, S}(k, params.ℓ, params.d, n, N, kskpar, gswpar, ffter, monomial, btk)
end

abstract type MKscheme{T, S} <: TFHEscheme{T, S} end

struct CCS{T, S} <: MKscheme{T, S}
    k::Int64
    n::Int64
    N::Int64
    a::TransCRS{S}
    kskpar::LEVparams_digit{T}
    unipar::Uniparams_digit{T}
    ffter::FFTransformer{S}
    monomial::Vector{TransNativePoly{S}}
    btk::Vector{BootKey_CCS{T, S}}
end

keygen_params(params::CCSparams{T, R, S}, ffter::FFTransformer{R}) where {T, R, S} =
    binary_lwekey(params.n, T), binary_ringkey(1, params.N, T, ffter)

"""
party_keygen outputs LWE key, RLWE key, and party-wise bootstrapping key.
"""
function party_keygen(a::CRS, params::CCSparams{T, R, S}) where {T, R, S}
    N = params.N

    kskpar = LEVparams_digit{T}(params.f, params.logD)
    unipar = Uniparams_digit{T}(1, params.l_uni, params.logB_uni)

    ffter_keygen = FFTransformer{R}(N, bits(T))
    ffter = FFTransformer{S}(N, bits(T))

    lwekey, ringkey = keygen_params(params, ffter_keygen)

    lwekey, ringkey, BootKey_CCS(lwekey, ringkey, kskpar, params.α, a, unipar, params.β, ffter_keygen, ffter)
end

"""
setup outputs Scheme.
"""
function setup(a::CRS, btk::Vector{BootKey_CCS{T, S}}, params::CCSparams{T, R, S}) where {T, R, S}
    k, n, N = params.k, params.n, params.N
    ffter = FFTransformer{S}(N, bits(T))
    kskpar = LEVparams_digit{T}(params.f, params.logD)
    unipar = Uniparams_digit{T}(1, params.l_uni, params.logB_uni)
    monomial = getmonomial(T, ffter)

    CCS{T, S}(k, n, N, fft.(a, Ref(ffter)), kskpar, unipar, ffter, monomial, btk)
end

abstract type KMSScheme{T, R, S} <: MKscheme{R, S} where T end

struct KMS{T, R, S} <: KMSScheme{T, R, S}
    k::Int64
    n::Int64
    N::Int64
    a::TransCRS{S}
    kskpar::LEVparams_digit{T}
    ffter::FFTransformer{S}
    monomial::Vector{TransNativePoly{S}}
    btk::Vector{BootKey_KMS{T, R, S}}
end

keygen_params(params::KMSparams{T, R, U, S}, ffter::FFTransformer{U}) where {T, R, U, S} =
    binary_lwekey(params.n, T), binary_ringkey(1, params.N, R, ffter), binary_ringkey(1, params.N, R, ffter)

"""
party_keygen outputs LWE key, RLWE key, and party-wise bootstrapping key.
"""
function party_keygen(a::CRS, params::KMSparams{T, R, U, S}) where {T, R, U, S}
    N = params.N

    kskpar = LEVparams_digit{T}(params.f, params.logD)
    gswpar = GSWparams_digit{R}(1, params.l_gsw, params.logB_gsw)
    levpar = LEVparams_digit{R}(params.l_lev, params.logB_lev)
    unipar = Uniparams_digit{R}(1, params.l_uni, params.logB_uni)

    ffter_keygen = FFTransformer{U}(N, bits(R))
    ffter = FFTransformer{S}(N, bits(R))

    lwekey, gswkey, unikey = keygen_params(params, ffter_keygen)

    lwekey, gswkey, unikey, BootKey_KMS(lwekey, gswkey, unikey, kskpar, params.α, levpar, gswpar, params.β, a, unipar, ffter_keygen, ffter)
end

"""
setup outputs Scheme.
"""
function setup(a::CRS, btk::Vector{BootKey_KMS{T, R, S}}, params::KMSparams{T, R, U, S}) where {T, R, U, S}
    k, n, N = params.k, params.n, params.N
    ffter = FFTransformer{S}(N, bits(R))
    kskpar = LEVparams_digit{T}(params.f, params.logD)
    monomial = getmonomial(T, ffter)

    KMS{T, R, S}(k, n, N, fft(a, ffter), kskpar, ffter, monomial, btk)
end

struct KMS_block{T, R, S} <: KMSScheme{T, R, S}
    k::Int64
    ℓ::Int64
    d::Int64
    n::Int64
    N::Int64
    a::TransCRS{S}
    kskpar::LEVparams_digit{T}
    ffter::FFTransformer{S}
    monomial::Vector{TransNativePoly{S}}
    btk::Vector{BootKey_KMS_block{T, R, S}}
end

keygen_params(params::KMSparams_block{T, R, U, S}, ffter::FFTransformer{U}) where {T, R, U, S} = begin
    lwekey = block_binary_lwekey(params.d, params.ℓ, T)
    gswkey = binary_ringkey(1, params.N, R, ffter)
    unikey = partial_ringkey(1, params.N, R, lwekey, ffter)
    lwekey, gswkey, unikey
end

"""
party_keygen outputs LWE key, RLWE key, and party-wise bootstrapping key.
"""
function party_keygen(a::CRS, params::KMSparams_block{T, R, U, S}) where {T, R, U, S}
    N = params.N

    kskpar = LEVparams_digit{T}(params.f, params.logD)
    gswpar = GSWparams_digit{R}(1, params.l_gsw, params.logB_gsw)
    levpar = LEVparams_digit{R}(params.l_lev, params.logB_lev)
    unipar = Uniparams_digit{R}(1, params.l_uni, params.logB_uni)

    ffter_keygen = FFTransformer{U}(N, bits(R))
    ffter = FFTransformer{S}(N, bits(R))

    lwekey, gswkey, unikey = keygen_params(params, ffter_keygen)

    lwekey, gswkey, unikey, BootKey_KMS_block(lwekey, gswkey, unikey, kskpar, params.α, levpar, gswpar, params.β, a, unipar, ffter_keygen, ffter, params.ℓ, params.d)
end

"""
setup outputs Scheme.
"""
function setup(a::CRS, btk::Vector{BootKey_KMS_block{T, R, S}}, params::KMSparams_block{T, R, U, S}) where {T, R, U, S}
    k, ℓ, d, N = params.k, params.ℓ, params.d, params.N
    ffter = FFTransformer{S}(N, bits(R))
    kskpar = LEVparams_digit{T}(params.f, params.logD)
    monomial = getmonomial(T, ffter)

    KMS_block{T, R, S}(k, ℓ, d, ℓ * d, N, fft(a, ffter), kskpar, ffter, monomial, btk)
end

function lwe_encrypt(m::Integer, key::LWEkey{T}, params::TFHEparams_bin) where T
    n = params.n
    b = unsigned(round(signed(T), gaussian(params.α)))
    a = rand(ChaCha20Stream(), T, n)
    μ = T(2) * T(m) - T(1)
    b += T(-reduce(+, a .* key.key) + μ << (bits(T) - 3))
    LWE(b, a)
end

function lwe_encrypt(m::Integer, key::LWEkey{T}, params::TFHEparams_block) where T
    n = params.ℓ * params.d
    b = unsigned(round(signed(T), gaussian(params.α)))
    a = rand(ChaCha20Stream(), T, n)
    μ = T(2) * T(m) - T(1)
    b += T(-reduce(+, a .* key.key) + μ << (bits(T) - 3))
    LWE(b, a)
end

function lwe_ith_encrypt(m::Integer, i::Int64, key::LWEkey{T}, params::KMSparams_block) where T
    n = params.ℓ * params.d
    b = unsigned(round(signed(T), gaussian(params.α)))
    a = vcat(zeros(T, n * (i-1)), rand(ChaCha20Stream(), T, n), zeros(T, n * (params.k-i)))
    μ = T(2) * T(m) - T(1)
    b += T(-reduce(+, a[n*(i-1)+1:n*i] .* key.key) + μ << (bits(T) - 3))
    LWE(b, a)
end

function lwe_ith_encrypt(m::Integer, i::Int64, key::LWEkey{T}, params::MKTFHEparams) where T
    n = params.n
    b = unsigned(round(signed(T), gaussian(params.α)))
    a = vcat(zeros(T, n * (i-1)), rand(ChaCha20Stream(), T, n), zeros(T, n * (params.k-i)))
    μ = T(2) * T(m) - T(1)
    b += T(-reduce(+, a[n*(i-1)+1:n*i] .* key.key) + μ << (bits(T) - 3))
    LWE(b, a)
end

lwe_decrypt(lwe::LWE{T}, key::LWEkey{T}) where T = 
    divbits(lwe.b + reduce(+, key.key .* lwe.a), bits(T) - 3) == 1

function lwe_decrypt(lwe::LWE{T}, keys::Vector{LWEkey{T}}, params::KMSparams_block) where T
    b, n = lwe.b, params.ℓ * params.d
    for i = eachindex(keys)
        b += reduce(+, keys[i].key .* lwe.a[(i-1)*n+1:i*n])
    end

    b < (T(1) << (bits(T) - 1))
end
    
function lwe_decrypt(lwe::LWE{T}, keys::Vector{LWEkey{T}}, params::MKTFHEparams) where T
    b, n = lwe.b, params.n
    for i = eachindex(keys)
        b += reduce(+, keys[i].key .* lwe.a[(i-1)*n+1:i*n])
    end

    b < (T(1) << (bits(T) - 1))
end

CRS(params::MKTFHEparams{T, R, S}) where {T, R, S} =
    [randnativepoly(params.N, T) for _ = 1 : params.l_uni]