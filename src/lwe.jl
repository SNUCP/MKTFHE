struct LweParams
    size :: Int
end


struct LweKey
    params :: LweParams
    key :: Array{Int32, 1}

    function LweKey(rng::AbstractRNG, params::LweParams)
        new(params, rand_uniform_bool(rng, params.size))
    end

    LweKey(params::LweParams, key::Array{Int32, 1}) = new(params, key)
end


"""
A structure representing an encrypted bit.
"""
mutable struct LweSample
    params :: LweParams
    a :: Array{Torus32, 1} # the n coefs of the mask
    b :: Torus32
    current_variance :: Float64 # average noise of the sample

    LweSample(params::LweParams, a::Array{Torus32, 1}, b::Torus32, cv::Float64) =
        new(params, a, b, cv)
end


Base.Broadcast.broadcastable(lwe::LweSample) = (lwe,)


"""
This function encrypts message by using key, with stdev alpha
"""
function lwe_encrypt(rng::AbstractRNG, message::Torus32, alpha::Float64, key::LweKey)
    lwe_params = key.params
    a = rand_uniform_torus32(rng, lwe_params.size)
    b = rand_gaussian_torus32(rng, message, alpha) + reduce(+, a .* key.key)
    LweSample(lwe_params, a, b, alpha^2)
end


"""
This function encrypts a message by using key and a given noise value
"""
function lwe_encrypt(
        rng::AbstractRNG, message::Torus32, noise::Float64, alpha::Float64, key::LweKey)
    lwe_params = key.params
    a = rand_uniform_torus32(rng, lwe_params.size)
    b = message + dtot32(noise) + reduce(+, a .* key.key)
    LweSample(lwe_params, a, b, alpha^2)
end


# This function computes the phase of sample by using key : phi = b - a.s
lwe_phase(x::LweSample, key::LweKey) = x.b - reduce(+, x.a .* key.key)


# result = (0,mu)
lwe_noiseless_trivial(mu::Torus32, params::LweParams) =
    LweSample(params, zeros(Torus32, params.size), mu, 0.)


Base.:+(x::LweSample, y::LweSample) =
    LweSample(x.params, x.a .+ y.a, x.b + y.b, x.current_variance + y.current_variance)


Base.:-(x::LweSample, y::LweSample) =
    LweSample(x.params, x.a .- y.a, x.b - y.b, x.current_variance + y.current_variance)

Base.:-(x::LweSample) = LweSample(x.params, -x.a, -x.b, x.current_variance)


function Base.:*(x::LweSample, y::Integer)
    ty = Torus32(y) # to make multiplication preserve the type
    LweSample(x.params, x.a .* ty, x.b * ty, x.current_variance * y^2)
end

Base.:*(x::Integer, y::LweSample) = y * x
