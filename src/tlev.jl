mutable struct TLevSample
    tlev_params :: TGswParams
    tlwe_params :: TLweParams
    samples :: Array{TLweSample, 1} # array of size (decomp_length, mask_size+1)

    TLevSample(tlev_params::TGswParams, samples::Array{TLweSample, 1}) =
        new(tlev_params, samples[1].params, samples)
end


Base.:+(x::TLevSample, y::TLevSample) =
    TLevSample(x.tlev_params, [x.samples[i] + y.samples[i] for i = 1 : x.tlev_params.decomp_length])


Base.:-(x::TLevSample, y::TLevSample) =
    TLevSample(x.tlev_params, [x.samples[i] - y.samples[i] for i = 1 : x.tlev_params.decomp_length])


TLevAddTo!(x::TLevSample, y::TLevSample) =
    x.samples .+= y.samples


TLevSubTo!(x::TLevSample, y::TLevSample) =
    x.samples .-= y.samples


struct TransformedTLevSample
    tlev_params :: TGswParams
    tlwe_params :: TLweParams
    samples :: Array{TransformedTLweSample, 1} # array of size (decomp_length, mask_size+1)

    TransformedTLevSample(tlev_params::TGswParams, samples::Array{TransformedTLweSample, 1}) =
        new(tlev_params, samples[1].params, samples)
end


function tlev_trivial_int(levpar::TGswParams, lwepar::TLweParams, message::Union{Int32, Int64})
    len = levpar.decomp_length
    deg = lwepar.polynomial_degree
    temp = lwepar.is32 ? 
        TLevSample(levpar, [tlwe_noiseless_trivial(zero_torus_polynomial(deg), lwepar) for i in 1:len]) :
        TLevSample(levpar, [tlwe_noiseless_trivial(zero_torus64_polynomial(deg), lwepar) for i in 1:len])
    tlev_add_gadget_times_message(temp, message)
end


function tlev_add_gadget_times_message(sample::TLevSample, message::Union{Int32, Int64})
    mask_size = sample.tlwe_params.mask_size
    decomp_length = sample.tlev_params.decomp_length
    gadget = sample.tlev_params.gadget_values

    samples = [
        TLweSample(
            sample.tlwe_params,
            sample.samples[i].a,
            sample.samples[i].current_variance)
        for i in 1:decomp_length]

    for i in 1:decomp_length
        samples[i].a[mask_size+1] += message*gadget[i]
    end

    TLevSample(sample.tlev_params, samples)
end


mul_by_monomial(x::TLevSample, shift::Integer) =
    TLevSample(x.tlev_params, mul_by_monomial.(x.samples, shift))


forward_transform(source::TLevSample) =
    TransformedTLevSample(source.tlev_params, forward_transform.(source.samples))


function tlev_extern_mul(c::TorusPolynomial, lev::TransformedTLevSample)
    dec_c = decompose(c, lev.tlev_params)
    tr_dec_c = forward_transform.(dec_c, lev.tlwe_params.is32)
    inverse_transform(sum(lev.samples .* tr_dec_c))
end


function tgsw_intern_mul(accum::TLevSample, gsw::TransformedTGswSample)
    samples = Array{TLweSample, 1}(undef, accum.tlev_params.decomp_length)
    for i in eachindex(samples)
        samples[i] = tgsw_extern_mul(accum.samples[i], gsw)
    end

    TLevSample(accum.tlev_params, samples)
end
