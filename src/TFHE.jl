module TFHE

using Random: AbstractRNG
using LinearAlgebra: mul!
using GenericFFT: plan_fft, plan_ifft, Plan
import DarkIntegers: mul_by_monomial
using DarkIntegers: Polynomial, negacyclic_modulus
using DoubleFloats: Double64

include("numeric-functions.jl")

include("polynomials.jl")

include("lwe.jl")

include("tlwe.jl")

include("tgsw.jl")

include("tlev.jl")

include("keyswitch.jl")

include("bootstrap.jl")

include("api.jl")
export make_key_pair
export LweSample
export SecretKey
export CloudKey
export encrypt
export decrypt
export tfhe_parameters_80
export tfhe_parameters_128

export SecretKey_new
export CloudKey_new

include("gates.jl")
export gate_nand
export gate_or
export gate_and
export gate_xor
export gate_xnor
export gate_not
export gate_constant
export gate_nor
export gate_andny
export gate_andyn
export gate_orny
export gate_oryn
export gate_mux

include("mk_internals.jl")

include("new_mk_internals.jl")
export test

include("mk_api.jl")
export SharedKey
export CloudKeyPart
export MKCloudKey
export mk_encrypt
export mk_decrypt
export mktfhe_parameters_2party
export mktfhe_parameters_4party
export mktfhe_parameters_8party
export mktfhe_parameters_16party

export SharedKey_new
export CloudKeyPart_new
export MKCloudKey_new
export mk_encrypt_new
export mk_decrypt_new
export mktfhe_parameters_2party_new
export mktfhe_parameters_2party_fast
export mktfhe_parameters_4party_new
export mktfhe_parameters_4party_fast
export mktfhe_parameters_8party_new
export mktfhe_parameters_8party_fast
export mktfhe_parameters_16party_new
export mktfhe_parameters_16party_fast
export mktfhe_parameters_32party_new
export mktfhe_parameters_32party_fast

include("mk_gates.jl")
export mk_gate_nand

include("new_mk_gates.jl")
export mk_gate_nand_new

end
