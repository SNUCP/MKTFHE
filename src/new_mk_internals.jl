struct BootstrapKeyPart_new

    uni_params :: TGswParams
    tlev_params :: TGswParams
    tgsw_params :: TGswParams
    tlwe_params :: TLweParams
    key_uni_enc :: MKTGswUESample
    gsw_key :: Array{TGswSample, 1}
    public_key :: PublicKey
    key_size :: Int

    function BootstrapKeyPart_new(
            rng, lwe_key::LweKey, tlwe_key::TLweKey, alpha_gsw::Float64, alpha_uni::Float64,
            shared_key::SharedKey, public_key::PublicKey,
            uni_params::TGswParams, tlev_params::TGswParams, tgsw_params::TGswParams)

        tlwe_params = shared_key.tlwe_params
        in_out_params = lwe_key.params

        rand_key = TLweKey(rng, tlwe_key.params)

        key_size = in_out_params.size
        gsw_enc = [
            tgsw_encrypt(rng, lwe_key.key[j], alpha_gsw, rand_key, tgsw_params)
            for j in 1:key_size]

        uni_enc = mk_tgsw_encrypt(rng, rand_key.key[1], alpha_uni, tlwe_key, shared_key, public_key)

        new(uni_params, tlev_params, tgsw_params, tlwe_params, uni_enc, gsw_enc, public_key, key_size)
    end
end


struct MKBootstrapKey_new

    uni_params :: TGswParams
    tlev_params :: TGswParams
    tgsw_params :: TGswParams
    tlwe_params :: TLweParams
    uni_key :: Array{MKTransformedTGswUESample, 1}
    gsw_key :: Array{TransformedTGswSample, 2}
    public_key :: Array{TransformedPublicKey}
    shared_key :: TransformedSharedKey

    function MKBootstrapKey_new(parts::Array{BootstrapKeyPart_new, 1}, sk::SharedKey)
        parties = length(parts)

        trans_uni = forward_transform.([parts[i].key_uni_enc for i in 1 : parties])
        trans_gsw = [forward_transform(parts[i].gsw_key[j]) for j in 1:parts[1].key_size, i in 1:parties]
        trans_pk = forward_transform.([parts[i].public_key for i in 1 : parties])
        trans_sk = forward_transform(sk)

        new(parts[1].uni_params, parts[1].tlev_params, parts[1].tgsw_params, parts[1].tlwe_params, trans_uni, trans_gsw, trans_pk, trans_sk)
    end
end


function UniProduct_new(
    acc::MKTLweSample, unienc::MKTransformedTGswUESample, pk::Array{TransformedPublicKey},
    sk::TransformedSharedKey, party::Int, parties::Int)

    tgsw_params = unienc.tgsw_params
    tlwe_params = acc.params
    declen = tgsw_params.decomp_length
    is32 = tlwe_params.is32

    dec_a = hcat(decompose.(acc.a, Ref(tgsw_params)))
    dec_b = decompose(acc.b, tgsw_params)

    tr_dec_a = is32 ? Array{TransformedTorusPolynomial}(undef, parties, declen) :
                      Array{TransformedTorusPolynomial64}(undef, parties, declen)
    tr_dec_b = is32 ? Array{TransformedTorusPolynomial}(undef, declen) :
                      Array{TransformedTorusPolynomial64}(undef, declen) 

    for i = 1 : parties
        tr_dec_a[i, :] = forward_transform.(dec_a[i], is32)
    end
    tr_dec_b = forward_transform.(dec_b, is32)

    u = Array{TorusPolynomial}(undef, parties)
    for i = 1 : parties
        u[i] = inverse_transform(sum([tr_dec_a[i, l] * unienc.d[l] for l in 1: declen]), is32)
    end
    u0 = inverse_transform(sum([tr_dec_b[l] * unienc.d[l] for l in 1:declen]), is32)

    v = inverse_transform(sum([tr_dec_a[i, l] * pk[i].b[l] for i in 1 : parties, l in 1: declen]), is32)
    v -= inverse_transform(sum([tr_dec_b[l] * sk.a[l] for l in 1: declen]), is32)

    dec_v = decompose(v, tgsw_params)
    tr_dec_v = forward_transform.(dec_v, is32)

    w0 = inverse_transform(sum([tr_dec_v[l] * unienc.f0[l] for l in 1:declen]), is32)
    w1 = inverse_transform(sum([tr_dec_v[l] * unienc.f1[l] for l in 1:declen]), is32)

    anew = u
    anew[party] += w1
    bnew = u0 + w0

    MKTLweSample(tlwe_params, anew, bnew, .0)
end


function mk_mux_rotate_new(
        accum::TLevSample, bki::TransformedTGswSample, barai::Int32)

    temp = mul_by_monomial(accum, barai) - accum
    accum + tgsw_intern_mul(temp, bki)
end


function mk_lev_tlwe_mul(accum::MKTLweSample, lev::TLevSample, rlk::MKTransformedTGswUESample, party::Int, 
                         parties::Int, pk::Array{TransformedPublicKey}, sk::TransformedSharedKey)
    N = accum.params.polynomial_degree
    zeropoly = accum.params.is32 ? zero_torus_polynomial : zero_torus64_polynomial

    e_a = [zeropoly(N) for i = 1 : parties]
    f_a = [zeropoly(N) for i = 1 : parties]
    fftlev = forward_transform(lev)

    for i = 1 : party-1
        temp = tlev_extern_mul(accum.a[i], fftlev)
        e_a[i] = temp.a[1]
        f_a[i] = temp.a[2]
    end
    temp = tlev_extern_mul(accum.b, fftlev)
    e_b = temp.a[1]
    f_b = temp.a[2]

    e = MKTLweSample(accum.params, e_a, e_b, 0.)
    f = MKTLweSample(accum.params, f_a, f_b, 0.)

    f - UniProduct_new(e, rlk, pk, sk, party, parties)
end


function mk_ith_blind_rotate(levpar::TGswParams, lwepar::TLweParams, 
            gsw_key::Array{TransformedTGswSample, 1}, bara::Array{Int32, 1})
    
    lwepar.is32 ? one = Int32(1) : one = Int64(1)
    levacc = tlev_trivial_int(levpar, lwepar, one)
    for i in eachindex(bara)
        barai = bara[i]
        if barai != 0
            levacc = mk_mux_rotate_new(levacc, gsw_key[i], barai)
        end
    end

    levacc
end


function mk_single_blind_rotate(tv::TorusPolynomial, lwepar::TLweParams, 
    gsw_key::Array{TransformedTGswSample, 1}, bara::Array{Int32, 1})

    acc = tlwe_noiseless_trivial(tv, lwepar)
    for i in eachindex(bara)
        barai = bara[i]
        if barai != 0
            acc = mux_rotate(acc, gsw_key[i], barai)
        end
    end

    acc
end


function mk_blind_rotate_new(accum::MKTLweSample, bk::MKBootstrapKey_new, bara::Array{Int32, 2})
    _, parties = size(bara)

    levkey = Array{TLevSample, 1}(undef, parties)

    for i in eachindex(levkey)
        levkey[i] = mk_ith_blind_rotate(bk.tlev_params, bk.tlwe_params, bk.gsw_key[:, i], bara[:, i])
        accum = mk_lev_tlwe_mul(accum, levkey[i], bk.uni_key[i], i, parties, bk.public_key, bk.shared_key)
    end

    accum
end


function mk_blind_rotate_new_v2(tv::TorusPolynomial, bk::MKBootstrapKey_new, bara::Array{Int32, 2})
    _, parties = size(bara)
    lwepar = bk.tlwe_params

    acc_1st = mk_single_blind_rotate(tv, bk.tlwe_params, bk.gsw_key[:, 1], bara[:, 1])
    e = mk_tlwe_noiseless_trivial(acc_1st.a[1], lwepar, parties)
    f = mk_tlwe_noiseless_trivial(acc_1st.a[2], lwepar, parties)

    accum = f - UniProduct_new(e, bk.uni_key[1], bk.public_key, bk.shared_key, 1, parties)

    levkey = Array{TLevSample, 1}(undef, parties-1)
    for i = 2 : parties
        levkey[i-1] = mk_ith_blind_rotate(bk.tlev_params, bk.tlwe_params, bk.gsw_key[:, i], bara[:, i])
        accum = mk_lev_tlwe_mul(accum, levkey[i-1], bk.uni_key[i], i, parties, bk.public_key, bk.shared_key)
    end

    accum
end


function mk_blind_rotate_and_extract_new(
        v::TorusPolynomial, bk::MKBootstrapKey_new, barb::Int32, bara::Array{Int32, 2})
    parties = size(bara, 2)
    testvectbis = mul_by_monomial(v, -barb)
    acc = mk_tlwe_noiseless_trivial(testvectbis, bk.tlwe_params, parties)
    acc = mk_blind_rotate_new(acc, bk, bara)
    bk.tlwe_params.is32 ? mk_tlwe_extract_sample(acc) : mk_tlwe_extract_sample_64(acc)
end


# This version skips the tlev partial bootstrapping for the first party.
function mk_blind_rotate_and_extract_new_v2(
        v::TorusPolynomial, bk::MKBootstrapKey_new, barb::Int32, bara::Array{Int32, 2})
    testvectbis = mul_by_monomial(v, -barb)
    acc = mk_blind_rotate_new_v2(testvectbis, bk, bara)
    bk.tlwe_params.is32 ? mk_tlwe_extract_sample(acc) : mk_tlwe_extract_sample_64(acc)
end


function mk_tlwe_extract_sample_64(x::MKTLweSample)
    a = t64tot32.(hcat([reverse_polynomial(p).coeffs for p in x.a]...))
    b = t64tot32(x.b.coeffs[1])
    lwe_params = LweParams(x.params.polynomial_degree * x.params.mask_size)
    MKLweSample(lwe_params, a, b, 0.)
end


function mk_bootstrap_wo_keyswitch_new(bk::MKBootstrapKey_new, mu::Union{Torus32, Torus64}, x::MKLweSample, fast_boot::Bool)

    p_degree = bk.tlwe_params.polynomial_degree
    barb = decode_message(x.b, p_degree * 2)
    bara = decode_message.(x.a, p_degree * 2)

    #the initial testvec = [mu,mu,mu,...,mu]
    testvect = torus_polynomial(repeat([mu], p_degree))

    fast_boot ? mk_blind_rotate_and_extract_new_v2(testvect, bk, barb, bara) : mk_blind_rotate_and_extract_new(testvect, bk, barb, bara)
end


function mk_bootstrap_new(bk::MKBootstrapKey_new, ks::Array{KeyswitchKey, 1}, mu::Union{Torus32, Torus64}, x::MKLweSample, fast_boot::Bool)
    u = mk_bootstrap_wo_keyswitch_new(bk, mu, x, fast_boot)
    mk_keyswitch(ks, u)
end