struct TGswParams
    decomp_length :: Int # decomposition length
    log2_base :: Int # log2(decomposition base)

    gadget_values :: Array{<:Union{Torus32, Torus64}, 1} # powers of the decomposition base
    offset :: Union{Int32, Int64} # decomposition offset

    is32 :: Bool

    function TGswParams(decomp_length::Int, log2_base::Int, is32::Bool)

        decomp_range = 1:decomp_length

        if is32
            torus = Torus32
            uint = UInt32
            bits = 32
        else
            torus = Torus64
            uint = UInt64
            bits = 64
        end

        # Nonzero entries of the gadget matrix.
        # 1/(base^(i+1)) as a Torus32 (the values of the gadget matrix)
        gadget_values = torus(2) .^ (bits .- decomp_range .* log2_base)

        # A precalculated offset value for decomposition of torus elements.
        # offset = base/2 * Sum{j=1..decomp_length} 2^(32 - j * bs_log2_base)
        offset = signed(uint(sum(gadget_values) * (Int128(2)^(log2_base - 1))))

        new(decomp_length, log2_base, gadget_values, offset, is32)
    end
end


struct TGswSample
    tgsw_params :: TGswParams
    tlwe_params :: TLweParams
    samples :: Array{TLweSample, 2} # array of size (decomp_length, mask_size+1)

    TGswSample(tgsw_params::TGswParams, samples::Array{TLweSample, 2}) =
        new(tgsw_params, samples[1].params, samples)
end


struct TransformedTGswSample
    tgsw_params :: TGswParams
    tlwe_params :: TLweParams
    samples :: Array{TransformedTLweSample, 2} # array of size (decomp_length, mask_size+1)

    TransformedTGswSample(tgsw_params::TGswParams, samples::Array{TransformedTLweSample, 2}) =
        new(tgsw_params, samples[1].params, samples)
end


"""
Returns `sample + message * h`, where `h` is the gadget matrix
(`h_ijk = delta_jk / B^i`, `i=1...decomp_length`, `j,k=1...mask_size+1`,
`B` is the decomposition base).
The dimensions of the TLWE sample array in the TGSW sample correspond to the first two indices
(`i` and `j`), and the index `k` corresponds to the vector of polynomials in a single TLWE sample.
"""
function tgsw_add_gadget_times_message(sample::TGswSample, message::Union{Int32, Int64})
    mask_size = sample.tlwe_params.mask_size
    decomp_length = sample.tgsw_params.decomp_length
    gadget = sample.tgsw_params.gadget_values

    # Since the gadget matrix is block-diagonal, we avoid extra work by only doing the addition
    # where its elements are nonzero.
    # (A violation of the Liskov principle here, which could be avoided by using special
    # types for zero TLWE samples and zero polynomials, but that's just too much infrastructure
    # for one small function)
    samples = [
        TLweSample(
            sample.tlwe_params,
            [j == k ? sample.samples[i,j].a[k] + message * gadget[i] : sample.samples[i,j].a[k]
                for k in 1:mask_size+1],
            # TODO: (issue #7) calculate current_variance correctly
            sample.samples[i,j].current_variance)
        for i in 1:decomp_length, j in 1:mask_size+1]

    TGswSample(sample.tgsw_params, samples)
end


function tgsw_encrypt_zero(rng::AbstractRNG, alpha::Float64, key::TLweKey, params::TGswParams)
    mask_size = key.params.mask_size
    decomp_length = params.decomp_length
    TGswSample(
        params,
        [tlwe_encrypt_zero(rng, alpha, key) for i in 1:decomp_length, j in 1:mask_size+1])
end


function tgsw_encrypt(
        rng::AbstractRNG, message::Int32, alpha::Float64, key::TLweKey, params::TGswParams)
    tgsw_zero = tgsw_encrypt_zero(rng, alpha, key, params)
    tgsw_add_gadget_times_message(tgsw_zero, message)
end


"""
Given the decomposition length `l`, decompose each coefficient of the given polynomial
into `l` values and store them in `l` polynomials.

For each value `x` in the real torus, the decomposition procedure floors it
to a multiple of `1/B^l` (where `B == 2^log2_base` is the decomposition base)
and finds `l` values `a_i` in `[-B/2, B/2) such that `x = sum(a_i / B^i, i=1...l)`.
"""
function decompose(sample::TorusPolynomial, params::TGswParams)

    decomp_length = params.decomp_length
    log2_base = params.log2_base

    if params.is32
        int = Int32
        bit = 32
    else
        int = Int64
        bit = 64
    end

    mask = int((1 << params.log2_base) - 1)
    part_offset = int(1 << (params.log2_base - 1))
    offset = params.offset

    # Since we want results in the range `[-B/2, B/2)` instead of `[0, B)`,
    # we add an `offset = sum(base powers)`,
    # decompose normally by calculating remainders
    # (since our base powers are powers of 2, these are implemented as shifts and ANDs)
    # and subtract back part of the offset (`B/2`).

    [int_polynomial(
        @. (((sample.coeffs + offset) >> (bit - power * log2_base)) & mask) - part_offset)
        for power in 1:decomp_length]
end


forward_transform(source::TGswSample) =
    TransformedTGswSample(source.tgsw_params, forward_transform.(source.samples))


# External product (*): accum = gsw (*) accum
function tgsw_extern_mul(accum::TLweSample, gsw::TransformedTGswSample)
    dec_accum = hcat(decompose.(accum.a, Ref(gsw.tgsw_params))...)
    tr_dec_accum = forward_transform.(dec_accum, accum.params.is32)
    inverse_transform(sum(gsw.samples .* tr_dec_accum))
end
