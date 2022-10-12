"""
TFHE scheme parameters (single- or multi- party).
"""
struct SchemeParameters

    lwe_size :: Int
    lwe_noise_stddev :: Float64

    tlwe_polynomial_degree :: Int
    tlwe_mask_size :: Int
    tlwe_is32 :: Bool

    bs_decomp_length :: Int # bootstrap decomposition length
    bs_log2_base :: Int # bootstrap log2(decomposition_base)
    bs_noise_stddev :: Float64 # bootstrap standard deviation

    ks_decomp_length :: Int # keyswitch decomposition length
    ks_log2_base :: Int # keyswitch log2(decomposition base)
    ks_noise_stddev :: Float64 # keyswitch noise standard deviation

    max_parties :: Int
end


struct SchemeParameters_new

    lwe_size :: Int
    lwe_noise_stddev :: Float64

    tlwe_polynomial_degree :: Int
    tlwe_mask_size :: Int
    tlwe_is32 :: Bool

    gsw_decomp_length :: Int # gsw decomposition length
    gsw_log2_base :: Int # gsw log2(decomposition_base)
    gsw_noise_stddev :: Float64 # gsw standard deviation

    lev_decomp_length :: Int # lev decomposition length
    lev_log2_base :: Int # lev log2(decomposition_base)

    uni_decomp_length :: Int # unienc decomposition length
    uni_log2_base :: Int # unienc log2(decomposition_base)
    uni_noise_stddev :: Float64 # unienc standard deviation

    ks_decomp_length :: Int # keyswitch decomposition length
    ks_log2_base :: Int # keyswitch log2(decomposition base)
    ks_noise_stddev :: Float64 # keyswitch noise standard deviation

    max_parties :: Int
end


"""
    tfhe_parameters_80(; tlwe_mask_size::Int=1)

Creates a single-party [`SchemeParameters`](@ref) object to pass to [`SecretKey`](@ref),
with parameters set to provide ~80 bits of security.
"""
tfhe_parameters_80(; tlwe_mask_size::Int=1) = SchemeParameters(
    # Parameters from I. Chillotti, N. Gama, M. Georgieva, and M. Izabachene,
    # "Faster Fully Homomorphic Encryption: Bootstrapping in Less Than 0.1 Seconds"
    # In 2020, estimated to provide ~80 bits of security.

    # LWE parameters
    500,
    1/2^15 * sqrt(2 / pi),
    # TLWE parameters
    1024, tlwe_mask_size,
    # bootstrap parameters
    2, 10, 9e-9 * sqrt(2 / pi),
    # keyswitch parameters
    8, 2, 1/2^15 * sqrt(2 / pi),
    1 # Only used for single-party encryption
    )



"""
    tfhe_parameters_128(; tlwe_mask_size::Int=1)

Creates a single-party [`SchemeParameters`](@ref) object to pass to [`SecretKey`](@ref),
with parameters set to provide ~128 bits of security.
"""
tfhe_parameters_128(; tlwe_mask_size::Int=1) = SchemeParameters(
    # Parameters from CGGI2019.
    # In 2020, estimated to provide ~129 bits of security.

    # LWE parameters
    630,
    1/2^15,
    # TLWE parameters
    1024, tlwe_mask_size,
    # bootstrap parameters
    3, 7, 1/2^25,
    # keyswitch parameters
    8, 2, 1/2^15,
    1 # Only used for single-party encryption
    )


lwe_parameters(params::SchemeParameters) =
    LweParams(params.lwe_size)

tlwe_parameters(params::SchemeParameters) =
    TLweParams(params.tlwe_polynomial_degree, params.tlwe_mask_size, params.tlwe_is32)

tgsw_parameters(params::SchemeParameters) =
    TGswParams(params.bs_decomp_length, params.bs_log2_base, params.tlwe_is32)

keyswitch_parameters(params::SchemeParameters) =
    KeyswitchParameters(params.ks_decomp_length, params.ks_log2_base)

"""
For SchemeParameters_new
"""

lwe_parameters(params::SchemeParameters_new) =
    LweParams(params.lwe_size)

tlwe_parameters(params::SchemeParameters_new) =
    TLweParams(params.tlwe_polynomial_degree, params.tlwe_mask_size, params.tlwe_is32)

tgsw_parameters(params::SchemeParameters_new) =
    TGswParams(params.gsw_decomp_length, params.gsw_log2_base, params.tlwe_is32)

tlev_parameters(params::SchemeParameters_new) =
    TGswParams(params.lev_decomp_length, params.lev_log2_base, params.tlwe_is32)

uni_parameters(params::SchemeParameters_new) =
    TGswParams(params.uni_decomp_length, params.uni_log2_base, params.tlwe_is32)

keyswitch_parameters(params::SchemeParameters_new) =
    KeyswitchParameters(params.ks_decomp_length, params.ks_log2_base)


"""
    SecretKey(rng::AbstractRNG, params::SchemeParameters)

A TFHE secret key, used for encryption/decryption.
Currently the only official way to get an object to pass to `params`
is from [`tfhe_parameters`](@ref).
"""
struct SecretKey
    params :: SchemeParameters
    key :: LweKey

    function SecretKey(rng::AbstractRNG, params::SchemeParameters)
        lwe_key = LweKey(rng, lwe_parameters(params))
        new(params, lwe_key)
    end
end


struct SecretKey_new
    params :: SchemeParameters_new
    key :: LweKey

    function SecretKey_new(rng::AbstractRNG, params::SchemeParameters_new)
        lwe_key = LweKey(rng, lwe_parameters(params))
        new(params, lwe_key)
    end
end


Base.Broadcast.broadcastable(sk::SecretKey) = (sk,)


"""
    CloudKey(rng::AbstractRNG, secret_key::SecretKey)

A TFHE cloud (public) key, used for secure computations by a third party.
"""
struct CloudKey
    params :: SchemeParameters
    bootstrap_key :: BootstrapKey
    keyswitch_key :: KeyswitchKey

    function CloudKey(rng::AbstractRNG, secret_key::SecretKey)
        params = secret_key.params
        tlwe_key = TLweKey(rng, tlwe_parameters(params))

        bs_key = BootstrapKey(
            rng, params.bs_noise_stddev, secret_key.key, tlwe_key, tgsw_parameters(params))
        ks_key = KeyswitchKey(
            rng, params.ks_noise_stddev, keyswitch_parameters(params), secret_key.key, tlwe_key)

        new(secret_key.params, bs_key, ks_key)
    end
end


Base.Broadcast.broadcastable(ck::CloudKey) = (ck,)


"""
    make_key_pair(rng::AbstractRNG, params::Union{Nothing, SchemeParameters}=nothing)

Creates a pair of [`SecretKey`](@ref) and a corresponding [`CloudKey`](@ref).
If `params` is `nothing`, the default return value of [`tfhe_parameters`](@ref) is used.
"""
function make_key_pair(rng::AbstractRNG, params::Union{Nothing, SchemeParameters}=nothing)
    if params === nothing
        params = tfhe_parameters_80()
    end
    secret_key = SecretKey(rng, params)
    cloud_key = CloudKey(rng, secret_key)
    secret_key, cloud_key
end


"""
    encrypt(rng::AbstractRNG, key::SecretKey, message::Bool)

Encrypts a plaintext bit.
Returns a [`LweSample`](@ref) object.
"""
function encrypt(rng::AbstractRNG, key::SecretKey, message::Bool)
    alpha = key.params.lwe_noise_stddev
    lwe_encrypt(rng, encode_message(message ? 1 : -1, 8), alpha, key.key)
end


"""
    decrypt(key::SecretKey, sample::LweSample)

Decrypts an encrypted bit.
Returns a boolean.
"""
function decrypt(key::SecretKey, sample::LweSample)
    lwe_phase(sample, key.key) > 0
end
