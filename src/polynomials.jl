IntPolynomial = Polynomial{<:Union{Int32, Int64}, N} where N
TorusPolynomial = Polynomial{<:Union{Torus32, Torus64}, N} where N


int_polynomial(coeffs) = Polynomial(coeffs, negacyclic_modulus)
torus_polynomial(coeffs) = Polynomial(coeffs, negacyclic_modulus)


zero_torus_polynomial(len) = torus_polynomial(zeros(Torus32, len))
zero_torus64_polynomial(len) = torus_polynomial(zeros(Torus64, len))


struct TransformedTorusPolynomial
    coeffs :: Array{Complex{Float64}, 1}
end


Base.:+(x::TransformedTorusPolynomial, y::TransformedTorusPolynomial) =
    TransformedTorusPolynomial(x.coeffs .+ y.coeffs)


Base.:*(x::TransformedTorusPolynomial, y::TransformedTorusPolynomial) =
    TransformedTorusPolynomial(x.coeffs .* y.coeffs)

function PolyAddTo!(x::TransformedTorusPolynomial, y::TransformedTorusPolynomial)
    x.coeffs .+= y.coeffs
end

function PolyMulTo!(x::TransformedTorusPolynomial, y::TransformedTorusPolynomial)
    x.coeffs .*= y.coeffs
end


Base.broadcastable(x::TransformedTorusPolynomial) = Ref(x)


struct TransformedTorusPolynomial64
    # coeffs :: Array{Complex{Double64}, 1}
    coeffs :: Array{Complex{Float64}, 1}
end


Base.:+(x::TransformedTorusPolynomial64, y::TransformedTorusPolynomial64) =
    TransformedTorusPolynomial64(x.coeffs .+ y.coeffs)


Base.:*(x::TransformedTorusPolynomial64, y::TransformedTorusPolynomial64) =
    TransformedTorusPolynomial64(x.coeffs .* y.coeffs)

function PolyAddTo!(x::TransformedTorusPolynomial64, y::TransformedTorusPolynomial64)
    x.coeffs .+= y.coeffs
end

function PolyMulTo!(x::TransformedTorusPolynomial64, y::TransformedTorusPolynomial64)
    x.coeffs .*= y.coeffs
end


Base.broadcastable(x::TransformedTorusPolynomial64) = Ref(x)


"""
For the given p(x), calculates p(1/x) with the applcation of the corresponding modulus
(x^N - 1) for cyclic polynomials, (x^N+1) for negacyclic ones.
"""
function reverse_polynomial(p::Polynomial)
    new_coeffs = collect(reverse(p.coeffs))
    mul_by_monomial(Polynomial(new_coeffs, p.modulus), length(new_coeffs) + 1)
end


#=
Using tangent FFT for polynomial convolution (by multiplying them in the transformed space).
Sacrificing some readability for a big speed improvement.
=#


struct ForwardTransformPlan

    plan :: Plan
    coeffs :: Array{Complex{Float64}, 1}
    buffer :: Array{Complex{Float64}, 1}

    function ForwardTransformPlan(len::Int)
        @assert len % 2 == 0
        idx = collect(0:len÷2-1)
        coeffs = exp.((-2im * pi / len / 2) .* idx)
        buffer = Array{Complex{Float64}, 1}(undef, len ÷ 2)
        plan = plan_fft(buffer)
        new(plan, coeffs, buffer)
    end
end


struct ForwardTransformPlan64

    plan :: Plan
    # coeffs :: Array{Complex{Double64}, 1}
    # buffer :: Array{Complex{Double64}, 1}
    coeffs :: Array{Complex{Float64}, 1}
    buffer :: Array{Complex{Float64}, 1}

    function ForwardTransformPlan64(len::Int)
        @assert len % 2 == 0
        idx = collect(0:len÷2-1)
        coeffs = exp.((-2im * pi / len / 2) .* idx)
        # buffer = Array{Complex{Double64}, 1}(undef, len ÷ 2)
        buffer = Array{Complex{Float64}, 1}(undef, len ÷ 2)
        plan = plan_fft(buffer)
        new(plan, coeffs, buffer)
    end
end


struct InverseTransformPlan

    plan :: Plan
    coeffs :: Array{Complex{Float64}, 1}
    complex_buffer :: Array{Complex{Float64}, 1}
    int_buffer :: Array{Torus32, 1}

    function InverseTransformPlan(len::Int)
        @assert len % 2 == 0
        idx = collect(0:len÷2-1)
        coeffs = exp.((-2im * pi / len / 2) .* idx)
        complex_buffer = Array{Complex{Float64}, 1}(undef, len ÷ 2)
        int_buffer = Array{Torus32, 1}(undef, len)
        plan = plan_ifft(complex_buffer)
        new(plan, coeffs, complex_buffer, int_buffer)
    end
end


struct InverseTransformPlan64

    plan :: Plan
    # coeffs :: Array{Complex{Double64}, 1}
    # complex_buffer :: Array{Complex{Double64}, 1}
    coeffs :: Array{Complex{Float64}, 1}
    complex_buffer :: Array{Complex{Float64}, 1}
    int_buffer :: Array{Torus64, 1}

    function InverseTransformPlan64(len::Int)
        @assert len % 2 == 0
        idx = collect(0:len÷2-1)
        coeffs = exp.((-2im * pi / len / 2) .* idx)
        # complex_buffer = Array{Complex{Double64}, 1}(undef, len ÷ 2)
        complex_buffer = Array{Complex{Float64}, 1}(undef, len ÷ 2)
        int_buffer = Array{Torus64, 1}(undef, len)
        plan = plan_ifft(complex_buffer)
        new(plan, coeffs, complex_buffer, int_buffer)
    end
end


_forward_transform_plans = Dict{Int, ForwardTransformPlan}()
_inverse_transform_plans = Dict{Int, InverseTransformPlan}()

_forward_transform_plans64 = Dict{Int, ForwardTransformPlan64}()
_inverse_transform_plans64 = Dict{Int, InverseTransformPlan64}()


function get_forward_transform_plan(len::Int, is32::Bool)
    if is32
        if !haskey(_forward_transform_plans, len)
            p = ForwardTransformPlan(len)
            _forward_transform_plans[len] = p
            p
        else
            _forward_transform_plans[len]
        end
    else
        if !haskey(_forward_transform_plans64, len)
            p = ForwardTransformPlan64(len)
            _forward_transform_plans64[len] = p
            p
        else
            _forward_transform_plans64[len]
        end
    end
end


function get_inverse_transform_plan(len::Int, is32::Bool)
    if is32
        if !haskey(_inverse_transform_plans, len)
            p = InverseTransformPlan(len)
            _inverse_transform_plans[len] = p
            p
        else
            _inverse_transform_plans[len]
        end
    else
        if !haskey(_inverse_transform_plans64, len)
            p = InverseTransformPlan64(len)
            _inverse_transform_plans64[len] = p
            p
        else
            _inverse_transform_plans64[len]
        end
    end
end


function forward_transform(p::Union{IntPolynomial, TorusPolynomial}, is32::Bool)
    c = p.coeffs
    N = length(c)
    p = get_forward_transform_plan(N, is32)
    p.buffer .= (c[1:N÷2] .- im .* c[N÷2+1:end]) .* p.coeffs
    is32 ? TransformedTorusPolynomial(p.plan * p.buffer) : TransformedTorusPolynomial64(p.plan * p.buffer)
end


to_int32(x::Int64) = signed(trunc(UInt32, unsigned(x) & 0xffffffff))
to_int32(x::Float64) = to_int32(round(Int64, x))
to_int64(x::Int128) = signed(trunc(UInt64, unsigned(x) & 0xffffffffffffffff))
to_int64(x::Float64) = to_int64(round(Int128, x))
# to_int64(x::Double64) = to_int64(round(Int128, x))


function inverse_transform(x::Union{TransformedTorusPolynomial, TransformedTorusPolynomial64}, is32::Bool)

    c = x.coeffs
    len = length(c)
    N = length(c)*2
    p = get_inverse_transform_plan(N, is32)

    mul!(p.complex_buffer, p.plan, x.coeffs)
    p.complex_buffer .= conj.(p.complex_buffer) .* p.coeffs
    if is32
        p.int_buffer[1:len] .= to_int32.(real.(p.complex_buffer))
        p.int_buffer[len+1:end] .= to_int32.(imag.(p.complex_buffer))
    else
        p.int_buffer[1:len] .= to_int64.(real.(p.complex_buffer))
        p.int_buffer[len+1:end] .= to_int64.(imag.(p.complex_buffer))
    end

    torus_polynomial(copy(p.int_buffer))
end


function transformed_mul(x::IntPolynomial, y::TorusPolynomial, is32::Bool)
    inverse_transform(forward_transform(x, is32) * forward_transform(y, is32), is32)
end
