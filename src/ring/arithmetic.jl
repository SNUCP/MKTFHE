@inline native(x::Float64, mask::UInt32) = begin
    x -= floor(x * 2.3283064365386963e-10) * 4.294967296e9
    unsafe_trunc(UInt32, x)
end

@inline native(x::Float64, mask::UInt64) = begin
    x -= floor(x * 5.421010862427522e-20) * 1.8446744073709552e19
    unsafe_trunc(UInt64, x)
end

@inline native(x::MultiFloat{Float64, l}, mask::T) where {T<:Unsigned, l} = begin
    res = zero(T)
    @inbounds for i = 1 : l
        res += native(x._limbs[i], mask)
    end
    res
end

to_int(a::Unsigned) = signed(a)

@inline bits(T::Type) = T == UInt32 ? 32 : 64

@inline divbits(a::T, bit::Int64) where {T<:Unsigned} = begin
    maxbit = bits(T)
    carry = (a << (maxbit - bit)) >> (maxbit - 1)
    a >> bit + carry 
end
