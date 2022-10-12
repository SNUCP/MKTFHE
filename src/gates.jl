#=
Homomorphic gates.
Take LWE samples with message in {-1/8, 1/8}, noise<1/16,
return a (bootstrapped) LWE sample with message in {-1/8, 1/8}, noise<1/16
(that is, a negative phase encodes `false`, a positive phase encodes `true`).
=#


"""
    gate_nand(ck::CloudKey, x::LweSample, y::LweSample)

Applies the NAND gate to encrypted bits `x` and `y`.
Returns a [`LweSample`](@ref) object.
"""
function gate_nand(ck::CloudKey, x::LweSample, y::LweSample)
    result = lwe_noiseless_trivial(encode_message(1, 8), x.params) - x - y
    bootstrap(ck.bootstrap_key, ck.keyswitch_key, encode_message(1, 8), result)
end


"""
    gate_or(ck::CloudKey, x::LweSample, y::LweSample)

Applies the OR gate to encrypted bits `x` and `y`.
Returns a [`LweSample`](@ref) object.
"""
function gate_or(ck::CloudKey, x::LweSample, y::LweSample)
    result = lwe_noiseless_trivial(encode_message(1, 8), x.params) + x + y
    bootstrap(ck.bootstrap_key, ck.keyswitch_key, encode_message(1, 8), result)
end


"""
    gate_and(ck::CloudKey, x::LweSample, y::LweSample)

Applies the AND gate to encrypted bits `x` and `y`.
Returns a [`LweSample`](@ref) object.
"""
function gate_and(ck::CloudKey, x::LweSample, y::LweSample)
    result = lwe_noiseless_trivial(encode_message(-1, 8), x.params) + x + y
    bootstrap(ck.bootstrap_key, ck.keyswitch_key, encode_message(1, 8), result)
end


"""
    gate_xor(ck::CloudKey, x::LweSample, y::LweSample)

Applies the XOR gate to encrypted bits `x` and `y`.
Returns a [`LweSample`](@ref) object.
"""
function gate_xor(ck::CloudKey, x::LweSample, y::LweSample)
    result = lwe_noiseless_trivial(encode_message(1, 4), x.params) + (x + y) * 2
    bootstrap(ck.bootstrap_key, ck.keyswitch_key, encode_message(1, 8), result)
end


"""
    gate_xnor(ck::CloudKey, x::LweSample, y::LweSample)

Applies the XNOR gate to encrypted bits `x` and `y`.
Returns a [`LweSample`](@ref) object.
"""
function gate_xnor(ck::CloudKey, x::LweSample, y::LweSample)
    result = lwe_noiseless_trivial(encode_message(-1, 4), x.params) - (x + y) * 2
    bootstrap(ck.bootstrap_key, ck.keyswitch_key, encode_message(1, 8), result)
end


"""
    gate_not(ck::CloudKey, x::LweSample)

Applies the NOT gate to an encrypted bit `x`.
Returns a [`LweSample`](@ref) object.
Does not need to be bootstrapped.
"""
function gate_not(ck::CloudKey, x::LweSample)
    # Not bootstrapped, the cloud key is just for the sake of interface uniformity.
    -x
end


"""
    gate_constant(ck::CloudKey, value::Bool)

Returns a [`LweSample`](@ref) object representing the plaintext `value`.

!!! note

    The returned object is suitable to use with other gates, but it is *not* encrypted.
"""
function gate_constant(ck::CloudKey, value::Bool)
    lwe_noiseless_trivial(encode_message(value ? 1 : -1, 8), lwe_parameters(ck.params))
end


"""
    gate_nor(ck::CloudKey, x::LweSample, y::LweSample)

Applies the NOR gate to encrypted bits `x` and `y`.
Returns a [`LweSample`](@ref) object.
"""
function gate_nor(ck::CloudKey, x::LweSample, y::LweSample)
    result = lwe_noiseless_trivial(encode_message(-1, 8), x.params) - x - y
    bootstrap(ck.bootstrap_key, ck.keyswitch_key, encode_message(1, 8), result)
end


"""
    gate_andny(ck::CloudKey, x::LweSample, y::LweSample)

Applies the ANDNY (ANDNY(x, y) == AND(NOT(x), y)) gate to encrypted bits `x` and `y`.
Returns a [`LweSample`](@ref) object.
"""
function gate_andny(ck::CloudKey, x::LweSample, y::LweSample)
    result = lwe_noiseless_trivial(encode_message(-1, 8), x.params) - x + y
    bootstrap(ck.bootstrap_key, ck.keyswitch_key, encode_message(1, 8), result)
end


"""
    gate_andyn(ck::CloudKey, x::LweSample, y::LweSample)

Applies the ANDYN (ANDYN(x, y) == AND(x, NOT(y))) gate to encrypted bits `x` and `y`.
Returns a [`LweSample`](@ref) object.
"""
function gate_andyn(ck::CloudKey, x::LweSample, y::LweSample)
    result = lwe_noiseless_trivial(encode_message(-1, 8), x.params) + x - y
    bootstrap(ck.bootstrap_key, ck.keyswitch_key, encode_message(1, 8), result)
end


"""
    gate_orny(ck::CloudKey, x::LweSample, y::LweSample)

Applies the ORNY (ORNY(x, y) == OR(NOT(x), y)) gate to encrypted bits `x` and `y`.
Returns a [`LweSample`](@ref) object.
"""
function gate_orny(ck::CloudKey, x::LweSample, y::LweSample)
    result = lwe_noiseless_trivial(encode_message(1, 8), x.params) - x + y
    bootstrap(ck.bootstrap_key, ck.keyswitch_key, encode_message(1, 8), result)
end


"""
    gate_oryn(ck::CloudKey, x::LweSample, y::LweSample)

Applies the ORYN (ORYN(x, y) == OR(x, NOT(y))) gate to encrypted bits `x` and `y`.
Returns a [`LweSample`](@ref) object.
"""
function gate_oryn(ck::CloudKey, x::LweSample, y::LweSample)
    result = lwe_noiseless_trivial(encode_message(1, 8), x.params) + x - y
    bootstrap(ck.bootstrap_key, ck.keyswitch_key, encode_message(1, 8), result)
end


"""
    gate_mux(ck::CloudKey, x::LweSample, y::LweSample, z::LweSample)

Applies the MUX (MUX(x, y, z) == x ? y : z == OR(AND(x, y), AND(NOT(x), z)))
gate to encrypted bits `x`, `y` and `z`.
Returns a [`LweSample`](@ref) object.
"""
function gate_mux(ck::CloudKey, x::LweSample, y::LweSample, z::LweSample)

    # compute `AND(x, y)`
    t1 = lwe_noiseless_trivial(encode_message(-1, 8), x.params) + x + y
    u1 = bootstrap_wo_keyswitch(ck.bootstrap_key, encode_message(1, 8), t1)

    # compute `AND(NOT(x), z)`
    t2 = lwe_noiseless_trivial(encode_message(-1, 8), x.params) - x + z
    u2 = bootstrap_wo_keyswitch(ck.bootstrap_key, encode_message(1, 8), t2)

    # compute `OR(u1,u2)`
    t3 = lwe_noiseless_trivial(encode_message(1, 8), u1.params) + u1 + u2

    keyswitch(ck.keyswitch_key, t3)
end
