struct KeyswitchParameters
    decomp_length :: Int # decomposition length
    log2_base :: Int # log2(decomposition base)
end


struct KeyswitchKey

    params :: KeyswitchParameters
    in_lwe_params :: LweParams
    out_lwe_params :: LweParams
    key :: Array{LweSample, 3} # the keyswitch elements: a (base-1, n, l) matrix

    function KeyswitchKey(
            rng::AbstractRNG, alpha::Float64,
            params::KeyswitchParameters, out_key::LweKey, tlwe_key::TLweKey)

        in_key = extract_lwe_key(tlwe_key)
        in_lwe_params = in_key.params
        out_lwe_params = out_key.params

        lwe_size = in_lwe_params.size
        decomp_length = params.decomp_length
        log2_base = params.log2_base
        base = 1 << log2_base

        # Generate centred noises
        noise = rand_gaussian_float(rng, alpha, lwe_size, decomp_length, base - 1)
        noise .-= sum(noise) / length(noise) # recentre

        # generate the keyswitch key
        # ks[h,j,i] encodes k.s[i]/base^(j+1) where `s` is the secret key.
        # (We're not storing the values for `h == 0`
        # since they will not be used during keyswitching)
        message(i,j,h) = (in_key.key[i] * Int32(h)) << (32 - j * log2_base)
        ks = [
            lwe_encrypt(rng, message(i,j,h), noise[i,j,h], alpha, out_key)
            for h in 1:base-1, j in 1:decomp_length, i in 1:lwe_size]

        new(params, in_lwe_params, out_lwe_params, ks)
    end
end


function keyswitch(ks::KeyswitchKey, sample::LweSample)
    lwe_size = ks.in_lwe_params.size
    log2_base = ks.params.log2_base
    decomp_length = ks.params.decomp_length

    result = lwe_noiseless_trivial(sample.b, ks.out_lwe_params)

    # Round `a` to the closest multiple of `1/2^t`, where `t` is the precision.
    # Since the real torus elements are stored as integers, adding
    # `2^32 * 1/2^t / 2 == 2^(32 - log2_base * decomp_length - 1)` sets the highmost
    # `log2_base * decomp_length` bits to the desired state
    # (similar to adding 0.5 before dropping everything to the right of the decimal point
    # when rounding a floating-point number).
    prec_offset = one(Int32) << (32 - (1 + log2_base * decomp_length))
    aibar = sample.a .+ prec_offset

    # Binary decompose the higmost bits of each `a_i`
    # into `decomp_length` values of `log2_base` bits each.
    base = one(Int32) << log2_base
    mask = base - one(Int32)
    a = [
        (ai >> (32 - j * log2_base)) & mask
        for ai in aibar, j in 1:decomp_length]

    # Translate the message of the result sample by -sum(a[i].s[i])
    # where s is the secret embedded in keyswitch key.
    for i in 1:lwe_size
        for j in 1:decomp_length
            if a[i,j] != 0
                result -= ks.key[a[i,j],j,i]
            end
        end
    end

    result
end
