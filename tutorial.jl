include("./src/TFHE.jl")
using .TFHE
using Random


function int_to_bits(x::T) where T <: Integer
    [((x >> (i - 1)) & 1 != 0) for i in 1:(sizeof(T)*8)]
end


function bits_to_int(::Type{T}, x::Array{Bool, 1}) where T <: Integer
    result = zero(T)
    for i in 1:min(sizeof(T) * 8, length(x))
        result |= (x[i]<<(i-1))
    end
    result
end


function prepare()

    rng = MersenneTwister(123)
    secret_key, cloud_key = make_key_pair(rng)

    # generate encrypt the 16 bits of 2017
    plaintext1 = UInt16(2017)
    bits1 = int_to_bits(plaintext1)
    ciphertext1 = [encrypt(rng, secret_key, bits1[1])]

    # generate encrypt the 16 bits of 42
    plaintext2 = UInt16(42)
    bits2 = int_to_bits(plaintext2)
    ciphertext2 = [encrypt(rng, secret_key, bits2[1])]

    secret_key, cloud_key, ciphertext1, ciphertext2
end


# elementary full comparator gate that is used to compare the i-th bit:
#   input: ai and bi the i-th bit of a and b
#          lsb_carry: the result of the comparison on the lowest bits
#   algo: if (a==b) return lsb_carry else return b
function encrypted_compare_bit(ck::CloudKey, a::LweSample, b::LweSample, lsb_carry::LweSample)
    tmp = gate_xnor(cloud_key, a, b)
    gate_mux(ck, tmp, lsb_carry, a)
end

# this function compares two multibit words, and puts the max in result
function encrypted_minimum(ck::CloudKey, a::Array{LweSample, 1}, b::Array{LweSample, 1})

    nb_bits = length(a)

    # initialize the carry to 0
    tmps1 = gate_constant(ck, false)
    # run the elementary comparator gate n times
    for i in 1:nb_bits
        @time tmps1 = encrypted_compare_bit(ck, a[i], b[i], tmps1)
    end

    # tmps1 is the result of the comparaison: 0 if a is larger, 1 if b is larger
    # select the max and copy it to the result
    @time [gate_mux(ck, tmps1, b[i], a[i]) for i in 1:nb_bits]
end


function process(ck::CloudKey, a::Array{LweSample, 1}, b::Array{LweSample, 1})
    # do some operations on the ciphertexts: here, we will compute the
    # minimum of the two
    encrypted_minimum(ck, ciphertext1, ciphertext2)
end


function verify(secret_key, answer)
    # decrypt and rebuild the 16-bit plaintext answer
    bits = [decrypt(secret_key, answer[i]) for i in 1:length(answer)]
    int_answer = bits_to_int(UInt16, bits)

    println("Answer: $int_answer")
end


secret_key, cloud_key, ciphertext1, ciphertext2 = prepare()
answer = process(cloud_key, ciphertext1, ciphertext2)
verify(secret_key, answer)
