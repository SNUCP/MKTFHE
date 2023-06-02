abstract type RingPoly{T<:Unsigned} <: AbstractVector{T} end
abstract type TransPoly{T} <: AbstractVector{T} end

mutable struct NativePoly{T} <: RingPoly{T}
    coeffs::Vector{T}
    const N::Int64
end

add(x::NativePoly, y::NativePoly) = begin
    @assert x.N == y.N "Two polynomials should have the same length."
    NativePoly(x.coeffs + y.coeffs, x.N)
end

add!(x::NativePoly, y::NativePoly) = begin
    addto!(x, x, y)
end

addto!(res::NativePoly, x::NativePoly, y::NativePoly) = begin
    @assert x.N == y.N "Two polynomials should have the same length."
    @assert res.N == x.N "The output and input polynomials should have the same length."
    @. res.coeffs = x.coeffs + y.coeffs
end

sub(x::NativePoly, y::NativePoly) = begin
    @assert x.N == y.N "Two polynomials should have the same length."
    NativePoly((@. x.coeffs - y.coeffs), x.N)
end

sub!(x::NativePoly, y::NativePoly) = begin
    subto!(x, x, y)
end

subto!(res::NativePoly, x::NativePoly, y::NativePoly) = begin
    @assert x.N == y.N "Two polynomials should have the same length."
    @assert res.N == x.N "The output and input polynomials should have the same length."
    @. res.coeffs = x.coeffs - y.coeffs
end

Base.:*(x::T, p::NativePoly{T}) where {T<:Unsigned} =
    NativePoly(x * p.coeffs, p.N)

initialise!(p::NativePoly{T}) where {T<:Unsigned} = 
    @. p.coeffs = zero(T)

zeronativepoly(N::Integer, T::Type) = 
    NativePoly(zeros(T, N), N)

randnativepoly(N::Integer, T::Type) =
    NativePoly(rand(T, N), N)

buffnativepoly(N::Integer, T::Type) =
    NativePoly(Vector{T}(undef, N), N)

# Changing TransNativePoly to immutable struct gives 5~7% of performance raise, with 3x allocations.
# Julia allows mutating a Vector inside the immutable instances.
mutable struct TransNativePoly{T<:AbstractFloat} <: TransPoly{T}
    coeffs::Vector{Complex{T}}
    const N::Int64
end

add(x::TransNativePoly, y::TransNativePoly) = begin
    @assert x.N == y.N "Two polynomials should have the same length."
    TransNativePoly((@. x.coeffs + y.coeffs), x.N)
end

add!(x::TransNativePoly, y::TransNativePoly) = begin
    addto!(x, x, y)
end

addto!(res::TransNativePoly, x::TransNativePoly, y::TransNativePoly) = begin
    @assert x.N == y.N "Two polynomials should have the same length."
    @assert res.N == x.N "The output and input polynomials should have the same length."
    @. res.coeffs = x.coeffs + y.coeffs
end

sub(x::TransNativePoly, y::TransNativePoly) = begin
    @assert x.N == y.N "Two polynomials should have the same length."
    TransNativePoly((@. x.coeffs - y.coeffs), x.N)
end

sub!(x::TransNativePoly, y::TransNativePoly) = begin
    subto!(x, x, y)
end

subto!(res::TransNativePoly, x::TransNativePoly, y::TransNativePoly) = begin
    @assert x.N == y.N "Two polynomials should have the same length."
    @assert res.N == x.N "The output and input polynomials should have the same length."
    @. res.coeffs = x.coeffs - y.coeffs
end

mul(x::TransNativePoly, y::TransNativePoly) = begin
    @assert x.N == y.N "Two Polynomials should have the same length."
    TransNativePoly((@. x.coeffs * y.coeffs), x.N)
end

mul!(x::TransNativePoly, y::TransNativePoly) = begin
    multo!(x, x, y)
end

multo!(res::TransNativePoly, x::TransNativePoly, y::TransNativePoly) = begin
    @assert x.N == y.N "Two Polynomials should have the same length."
    @. res.coeffs = x.coeffs * y.coeffs
end

# This function is the bottleneck.
muladdto!(res::TransNativePoly, x::TransNativePoly, y::TransNativePoly) = begin
    @assert x.N == y.N "Two Polynomials should have the same length."
    @. res.coeffs += x.coeffs * y.coeffs
end

mulsubto!(res::TransNativePoly, x::TransNativePoly, y::TransNativePoly) = begin
    @assert x.N == y.N "Two Polynomials should have the same length."
    @. res.coeffs -= x.coeffs * y.coeffs
end

initialise!(p::TransNativePoly{T}) where {T<:AbstractFloat} = 
    @. p.coeffs = 0

zerotransnativepoly(N::Integer, T::Type) = 
    TransNativePoly(zeros(Complex{T}, N÷2), N)

randtransnativepoly(N::Integer, T::Type) =
    TransNativePoly(rand(Complex{T}, N÷2), N)

bufftransnativepoly(N::Integer, T::Type) =
    TransNativePoly(Vector{Complex{T}}(undef, N÷2), N)