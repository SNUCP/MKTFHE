struct TLweParams
    polynomial_degree :: Int # a power of 2: degree of the polynomials
    mask_size :: Int # number of polynomials in the mask

    is32 :: Bool

    function TLweParams(polynomial_degree::Int, mask_size::Int, is32::Bool)
        new(polynomial_degree, mask_size, is32)
    end
end


struct TLweKey
    params :: TLweParams
    key :: Array{IntPolynomial, 1} # the key (i.e k binary polynomials)

    function TLweKey(rng::AbstractRNG, params::TLweParams)
        key = params.is32 ?
            [int_polynomial(rand_uniform_bool(rng, params.polynomial_degree)) for i in 1:params.mask_size] : 
            [int_polynomial(rand_uniform_bool64(rng, params.polynomial_degree)) for i in 1:params.mask_size]
        new(params, key)
    end
end


# extractions Ring Lwe . Lwe
function extract_lwe_key(tlwe_key::TLweKey)
    tlwe_params = tlwe_key.params

    key = Int32.(vcat([poly.coeffs for poly in tlwe_key.key]...))

    LweKey(LweParams(tlwe_params.polynomial_degree * tlwe_params.mask_size), key)
end


struct TLweSample
    params :: TLweParams
    a :: Array{TorusPolynomial, 1} # array of length mask_size+1: mask + right term
    current_variance :: Float64 # avg variance of the sample

    TLweSample(params::TLweParams, a::Array{T, 1}, cv::Float64) where T <: TorusPolynomial =
        new(params, a, cv)
end


struct TransformedTLweSample
    params :: TLweParams
    a :: Array{<:Union{TransformedTorusPolynomial, TransformedTorusPolynomial64}, 1} # array of length mask_size+1: mask + right term
    current_variance :: Float64 # avg variance of the sample

    TransformedTLweSample(
            params::TLweParams, a::Array{<:Union{TransformedTorusPolynomial, TransformedTorusPolynomial64}, 1}, cv::Float64) =
        new(params, a, cv)
end


function tlwe_extract_sample(x::TLweSample)
    a = vcat([reverse_polynomial(p).coeffs for p in x.a[1:end-1]]...)
    b = x.a[end].coeffs[1]
    LweSample(LweParams(length(a)), a, b, 0.) # TODO: (issue #7) calculate the current variance
end


# create an homogeneous tlwe sample
function tlwe_encrypt_zero(rng::AbstractRNG, alpha::Float64, key::TLweKey)
    params = key.params
    polynomial_degree = params.polynomial_degree
    mask_size = params.mask_size

    if params.is32
        rand_uni = rand_uniform_torus32
        rand_gauss = rand_gaussian_torus32
        zero = Int32(0)
    else
        rand_uni = rand_uniform_torus64
        rand_gauss = rand_gaussian_torus64
        zero = Int64(0)
    end

    a_part = [torus_polynomial(rand_uni(rng, polynomial_degree)) for i in 1:mask_size]
    a_last = (
        int_polynomial(rand_gauss(rng, zero, alpha, polynomial_degree))
        + sum(transformed_mul.(key.key, a_part, params.is32)))
    TLweSample(params, vcat(a_part, [a_last]), alpha^2)
end


# result = (0,mu)
function tlwe_noiseless_trivial(mu::TorusPolynomial, params::TLweParams)
    a_part = params.is32 ? 
        [zero_torus_polynomial(params.polynomial_degree) for i in 1:params.mask_size] : 
        [zero_torus64_polynomial(params.polynomial_degree) for i in 1:params.mask_size]
    a_last = deepcopy(mu)
    TLweSample(params, vcat(a_part, [a_last]), 0.)
end


Base.:+(x::TLweSample, y::TLweSample) =
    TLweSample(x.params, x.a .+ y.a, x.current_variance + y.current_variance)


Base.:-(x::TLweSample, y::TLweSample) =
    TLweSample(x.params, x.a .- y.a, x.current_variance + y.current_variance)


mul_by_monomial(x::TLweSample, shift::Integer) =
    TLweSample(x.params, mul_by_monomial.(x.a, shift), x.current_variance)


forward_transform(x::TLweSample) =
    TransformedTLweSample(x.params, forward_transform.(x.a, x.params.is32), x.current_variance)


inverse_transform(x::TransformedTLweSample) =
    TLweSample(x.params, inverse_transform.(x.a, x.params.is32), x.current_variance)


# TODO: (issue #7) how to compute the variance correctly?
Base.:+(x::TransformedTLweSample, y::TransformedTLweSample) =
    TransformedTLweSample(x.params, x.a .+ y.a, x.current_variance + y.current_variance)


# TODO: (issue #7) how to compute the variance correctly?
Base.:*(x::TransformedTLweSample, y::Union{TransformedTorusPolynomial, TransformedTorusPolynomial64}) =
    TransformedTLweSample(x.params, x.a .* y, x.current_variance)
