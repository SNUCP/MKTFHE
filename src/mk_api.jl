"""
Multi-key TFHE parameters for 2 parties (a [`SchemeParameters`](@ref) object).
"""
mktfhe_parameters_2party = SchemeParameters(
    560, 3.05e-5, # LWE parameters
    1024, 1, true, # TLWE parameters
    3, 9, 3.72e-9, # bootstrap parameters
    8, 2, 3.05e-5, # keyswitch parameters
    2
    )

mktfhe_parameters_2party_new = SchemeParameters_new(
    560, 3.05e-5, # LWE parameters
    2048, 1, false, # TLWE parameters
    3, 13, 4.63e-18, # gsw parameters
    2, 7, # lev parameters
    2, 13, 4.63e-18, # unienc parameters
    8, 2, 3.05e-5, # keyswitch parameters
    2
    )

mktfhe_parameters_2party_fast = SchemeParameters_new(
    560, 3.05e-5, # LWE parameters
    2048, 1, false, # TLWE parameters
    3, 13, 4.63e-18, # gsw parameters
    2, 7, # lev parameters
    3, 10, 4.63e-18, # unienc parameters
    8, 2, 3.05e-5, # keyswitch parameters
    2
    )

"""
Multi-key TFHE parameters for 4 parties (a [`SchemeParameters`](@ref) object).
"""
mktfhe_parameters_4party = SchemeParameters(
    560, 3.05e-5, # LWE parameters
    1024, 1, true,# TLWE parameters
    4, 8, 3.72e-9, # bootstrap parameters
    8, 2, 3.05e-5, # keyswitch parameters
    4
    )

mktfhe_parameters_4party_new = SchemeParameters_new(
    560, 3.05e-5, # LWE parameters
    2048, 1, false, # TLWE parameters
    5, 8, 4.63e-18, # gsw parameters #37
    2, 8, # lev parameters
    5, 8, 4.63e-18, # unienc parameters
    8, 2, 3.05e-5, # keyswitch parameters
    4
    )

mktfhe_parameters_4party_fast = SchemeParameters_new(
    560, 3.05e-5, # LWE parameters
    2048, 1, false, # TLWE parameters
    5, 8, 4.63e-18, # gsw parameters #37
    2, 8, # lev parameters
    7, 6, 4.63e-18, # unienc parameters
    8, 2, 3.05e-5, # keyswitch parameters
    4
    )


"""
Multi-key TFHE parameters for 8 parties (a [`SchemeParameters`](@ref) object).
"""
mktfhe_parameters_8party = SchemeParameters(
    560, 3.05e-5, # LWE parameters
    1024, 1, true,# TLWE parameters
    5, 6, 3.72e-9, # bootstrap parameters
    8, 2, 3.05e-5, # keyswitch parameters
    8
    )


mktfhe_parameters_8party_new = SchemeParameters_new(
    560, 3.05e-5, # LWE parameters
    2048, 1, false, # TLWE parameters
    4, 11, 4.63e-18, # gsw parameters
    3, 6, # lev parameters
    8, 4, 4.63e-18, # unienc parameters
    8, 2, 3.05e-5, # keyswitch parameters
    8
    )

mktfhe_parameters_8party_fast = SchemeParameters_new(
    560, 3.05e-5, # LWE parameters
    2048, 1, false, # TLWE parameters
    4, 11, 4.63e-18, # gsw parameters
    3, 6, # lev parameters
    7, 4, 4.63e-18, # unienc parameters
    8, 2, 3.05e-5, # keyswitch parameters
    8
    )


"""
Multi-key TFHE parameters for 16 parties (a [`SchemeParameters`](@ref) object).
"""
mktfhe_parameters_16party = SchemeParameters(
    560, 3.05e-5, # LWE parameters
    1024, 1, true,# TLWE parameters
    12, 2, 3.72e-9, # bootstrap parameters
    8, 2, 3.05e-5, # keyswitch parameters
    16
    )


mktfhe_parameters_16party_new = SchemeParameters_new(
    560, 3.05e-5, # LWE parameters
    2048, 1, false, # TLWE parameters
    5, 9, 4.63e-18, # gsw parameters
    3, 6, # lev parameters
    9, 4, 4.63e-18, # unienc parameters
    8, 2, 3.05e-5, # keyswitch parameters
    16
    )

mktfhe_parameters_16party_fast = SchemeParameters_new(
    560, 3.05e-5, # LWE parameters
    2048, 1, false, # TLWE parameters
    5, 9, 4.63e-18, # gsw parameters
    3, 6, # lev parameters
    7, 4, 4.63e-18, # unienc parameters
    8, 2, 3.05e-5, # keyswitch parameters
    16
    )

"""
Multi-key TFHE parameters for 32 parties (a [`SchemeParameters`](@ref) object).
"""
mktfhe_parameters_32party_new = SchemeParameters_new(
    560, 3.05e-5, # LWE parameters
    2048, 1, false, # TLWE parameters
    6, 8, 4.63e-18, # gsw parameters
    3, 7, # lev parameters
    16, 2, 4.63e-18, # unienc parameters
    8, 2, 3.05e-5, # keyswitch parameters
    32
    )

mktfhe_parameters_32party_fast = SchemeParameters_new(
    560, 3.05e-5, # LWE parameters
    2048, 1, false, # TLWE parameters
    6, 8, 4.63e-18, # gsw parameters
    3, 7, # lev parameters
    16, 2, 4.63e-18, # unienc parameters
    8, 2, 3.05e-5, # keyswitch parameters
    32
    )

"""
    SharedKey(rng::AbstractRNG, params::SchemeParameters)

A shared key created by the server.
`params` is one of [`mktfhe_parameters_2party`](@ref), [`mktfhe_parameters_4party`](@ref),
[`mktfhe_parameters_8party`](@ref).
"""
function SharedKey(rng::AbstractRNG, params::SchemeParameters)
    # Resolving a circular dependency. SharedKey is used internally,
    # and we don't want to expose SchemeParameters there.
    tgsw_params = tgsw_parameters(params)
    tlwe_params = tlwe_parameters(params)
    SharedKey(rng, tgsw_params, tlwe_params)
end


function SharedKey_new(rng::AbstractRNG, params::SchemeParameters_new)
    # Resolving a circular dependency. SharedKey is used internally,
    # and we don't want to expose SchemeParameters there.
    uni_params = uni_parameters(params)
    tlwe_params = tlwe_parameters(params)
    SharedKey(rng, uni_params, tlwe_params)
end


"""
    CloudKeyPart(rng, secret_key::SecretKey, shared_key::SharedKey)

A part of the cloud (computation) key generated independently by each party
(since it involves their secret keys).
The `secret_key` is a [`SecretKey`](@ref) object created with the same parameter set
as the [`SharedKey`](@ref) object.
"""
struct CloudKeyPart
    params :: SchemeParameters
    bk_part :: BootstrapKeyPart
    ks :: KeyswitchKey

    function CloudKeyPart(rng, secret_key::SecretKey, shared_key::SharedKey)
        params = secret_key.params
        tgsw_params = tgsw_parameters(params)
        tlwe_key = TLweKey(rng, tlwe_parameters(params))
        pk = PublicKey(rng, tlwe_key, params.bs_noise_stddev, shared_key, tgsw_params)
        bk = BootstrapKeyPart(rng, secret_key.key, tlwe_key, params.bs_noise_stddev, shared_key, pk)
        ks = KeyswitchKey(
                rng, params.ks_noise_stddev, keyswitch_parameters(params),
                secret_key.key, tlwe_key)
        new(params, bk, ks)
    end
end


"""
    MKCloudKey(ck_parts::Array{CloudKeyPart, 1})

A full cloud key generated on the server out of parties' cloud key parts.
"""
struct MKCloudKey
    parties :: Int
    params :: SchemeParameters
    bootstrap_key :: MKBootstrapKey
    keyswitch_key :: Array{KeyswitchKey, 1}

    function MKCloudKey(ck_parts::Array{CloudKeyPart, 1}, shared_key::SharedKey)
        params = ck_parts[1].params
        parties = length(ck_parts)
        @assert parties <= params.max_parties

        bk = MKBootstrapKey([part.bk_part for part in ck_parts], shared_key)
        ks = [part.ks for part in ck_parts]

        new(parties, params, bk, ks)
    end
end

"""
New.
"""

struct CloudKeyPart_new
    params :: SchemeParameters_new
    bk_part :: BootstrapKeyPart_new
    ks :: KeyswitchKey

    function CloudKeyPart_new(rng, secret_key::SecretKey_new, shared_key::SharedKey)
        params = secret_key.params
        uni_params = uni_parameters(params)
        tlev_params = tlev_parameters(params)
        tgsw_params = tgsw_parameters(params)

        tlwe_key = TLweKey(rng, tlwe_parameters(params))

        pk = PublicKey(rng, tlwe_key, params.uni_noise_stddev, shared_key, uni_params)
        bk = BootstrapKeyPart_new(rng, secret_key.key, tlwe_key, params.gsw_noise_stddev, params.uni_noise_stddev,
                shared_key, pk, uni_params, tlev_params, tgsw_params)
        ks = KeyswitchKey(
                rng, params.ks_noise_stddev, keyswitch_parameters(params),
                secret_key.key, tlwe_key)

        new(params, bk, ks)
    end
end


struct MKCloudKey_new
    parties :: Int
    params :: SchemeParameters_new
    bootstrap_key :: MKBootstrapKey_new
    keyswitch_key :: Array{KeyswitchKey, 1}

    function MKCloudKey_new(ck_parts::Array{CloudKeyPart_new, 1}, shared_key::SharedKey)
        params = ck_parts[1].params
        parties = length(ck_parts)
        @assert parties <= params.max_parties

        bk = MKBootstrapKey_new([part.bk_part for part in ck_parts], shared_key)
        ks = [part.ks for part in ck_parts]

        new(parties, params, bk, ks)
    end
end


"""
    mk_encrypt(rng, secret_keys::Array{SecretKey, 1}, message::Bool)

Encrypts a plaintext bit using parties' secret keys.
Returns a [`MKLweSample`](@ref) object.
"""
function mk_encrypt(rng, secret_keys::Array{SecretKey, 1}, message::Bool)

    # TODO: (issue #6) encrypt separately for each party and share <a_i, s_i>?

    mu = encode_message(message ? 1 : -1, 8)

    params = secret_keys[1].params
    lwe_params = lwe_parameters(params)
    alpha = params.lwe_noise_stddev
    parties = length(secret_keys)

    a = hcat([rand_uniform_torus32(rng, lwe_params.size) for i in 1:parties]...)
    b = (rand_gaussian_torus32(rng, mu, alpha)
        + reduce(+, a .* hcat([secret_keys[i].key.key for i in 1:parties]...)))

    MKLweSample(lwe_params, a, b, alpha^2)
end


function mk_encrypt_new(rng, secret_keys::Array{SecretKey_new, 1}, message::Bool)

    # TODO: (issue #6) encrypt separately for each party and share <a_i, s_i>?

    mu = encode_message(message ? 1 : -1, 8)

    params = secret_keys[1].params
    lwe_params = lwe_parameters(params)
    alpha = params.lwe_noise_stddev
    parties = length(secret_keys)

    a = hcat([rand_uniform_torus32(rng, lwe_params.size) for i in 1:parties]...)
    b = (rand_gaussian_torus32(rng, mu, alpha)
        + reduce(+, a .* hcat([secret_keys[i].key.key for i in 1:parties]...)))

    MKLweSample(lwe_params, a, b, alpha^2)
end


"""
    mk_decrypt(secret_keys::Array{SecretKey, 1}, sample::MKLweSample)

Decrypts an encrypted bit using parties' secret keys.
Returns a boolean.
"""
function mk_decrypt(secret_keys::Array{SecretKey, 1}, sample::MKLweSample)
    # TODO: (issue #6) decrypt separately at each party and join phases?
    mk_lwe_phase(sample, [sk.key for sk in secret_keys]) > 0
end

function mk_decrypt_new(secret_keys::Array{SecretKey_new, 1}, sample::MKLweSample)
    # TODO: (issue #6) decrypt separately at each party and join phases?
    mk_lwe_phase(sample, [sk.key for sk in secret_keys]) > 0
end
