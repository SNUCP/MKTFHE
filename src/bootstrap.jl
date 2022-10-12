struct BootstrapKey
    tgsw_params :: TGswParams # params of the Gsw elems in bk. key: s"
    tlwe_params :: TLweParams
    key :: Array{TransformedTGswSample, 1} # the bootstrapping key (s->s")

    function BootstrapKey(
            rng::AbstractRNG, alpha::Float64, lwe_key::LweKey,
            tlwe_key::TLweKey, tgsw_params::TGswParams)

        tlwe_params = tlwe_key.params
        bk = [tgsw_encrypt(rng, key_elem, alpha, tlwe_key, tgsw_params) for key_elem in lwe_key.key]
        transformed_bk = forward_transform.(bk)

        new(tgsw_params, tlwe_params, transformed_bk)
    end
end


function mux_rotate(accum::TLweSample, bki::TransformedTGswSample, barai::Int32)
    # accum += BK_i * [(X^bar{a}_i-1) * accum]
    temp = mul_by_monomial(accum, barai) - accum
    accum + tgsw_extern_mul(temp, bki)
end


#=
 * multiply the accumulator by X^sum(bara_i.s_i)
 * @param accum the TLWE sample to multiply
 * @param bk An array of n TGSW FFT samples where bk_i encodes s_i
 * @param bara An array of n coefficients between 0 and 2N-1
=#
function blind_rotate(accum::TLweSample, bk::BootstrapKey, bara::Array{Int32, 1})
    for i in 1:length(bk.key)
        if bara[i] != 0
            accum = mux_rotate(accum, bk.key[i], bara[i])
        end
    end
    accum
end


#=
 * result = LWE(v_p) where p=barb-sum(bara_i.s_i) mod 2N
 * @param result the output LWE sample
 * @param v a 2N-elt anticyclic function (represented by a TorusPolynomial)
 * @param bk An array of n TGSW FFT samples where bk_i encodes s_i
 * @param barb A coefficients between 0 and 2N-1
 * @param bara An array of n coefficients between 0 and 2N-1
=#
function blind_rotate_and_extract(
        v::TorusPolynomial, bk::BootstrapKey, barb::Int32, bara::Array{Int32, 1})

    # testvector = X^{2N-barb}*v == X^{-barb}*v
    testvectbis = mul_by_monomial(v, -barb)

    accum = tlwe_noiseless_trivial(testvectbis, bk.tlwe_params)
    accum = blind_rotate(accum, bk, bara)
    tlwe_extract_sample(accum)
end


#=
 * result = LWE(mu) iff phase(x)>0, LWE(-mu) iff phase(x)<0
 * @param result The resulting LweSample
 * @param bk The bootstrapping + keyswitch key
 * @param mu The output message (if phase(x)>0)
 * @param x The input sample
=#
function bootstrap_wo_keyswitch(bk::BootstrapKey, mu::Torus32, x::LweSample)

    p_degree = bk.tlwe_params.polynomial_degree

    # Modulus switching
    bara = decode_message.(x.a, p_degree * 2)
    barb = decode_message(x.b, p_degree * 2)

    # the initial testvec = [mu,mu,mu,...,mu]
    testvect = torus_polynomial(repeat([mu], p_degree))

    # Bootstrapping rotation and extraction
    blind_rotate_and_extract(testvect, bk, barb, bara)
end


#=
 * result = LWE(mu) iff phase(x)>0, LWE(-mu) iff phase(x)<0
 * @param result The resulting LweSample
 * @param bk The bootstrapping + keyswitch key
 * @param mu The output message (if phase(x)>0)
 * @param x The input sample
=#
function bootstrap(bk::BootstrapKey, ks::KeyswitchKey, mu::Torus32, x::LweSample)
    u = bootstrap_wo_keyswitch(bk, mu, x)
    keyswitch(ks, u)
end
