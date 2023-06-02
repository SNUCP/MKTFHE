@inline native(x::Float64, mask::T) where {T<:Unsigned} = begin
    if x == .0
        zero(T)
    else
        shift = exponent(x) - 52
        if shift ≥ 0 
            x > 0 ? T(((reinterpret(UInt64, x) & 0x000fffffffffffff + 0x0010000000000000) << shift) & mask) :
                   -T(((reinterpret(UInt64, x) & 0x000fffffffffffff + 0x0010000000000000) << shift) & mask)
        else
            x > 0 ? T(((reinterpret(UInt64, x) & 0x000fffffffffffff + 0x0010000000000000) >> -shift) & mask) :
                   -T(((reinterpret(UInt64, x) & 0x000fffffffffffff + 0x0010000000000000) >> -shift) & mask)
        end
    end
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