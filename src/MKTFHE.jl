module MKTFHE

using Random: RandomDevice
using Base.Threads: @threads, @spawn, @sync
using MultiFloats
# MultiFloats package is the state-of-the-art multi-precision Julia package.

include("ring/arithmetic.jl")
include("ring/polynomial.jl")
include("ring/sampler.jl")
include("ring/fft.jl")

include("ciphertext/key.jl")
include("ciphertext/lwe.jl")
include("ciphertext/lev.jl")
include("ciphertext/gsw.jl")
include("ciphertext/unienc.jl")

include("tfhe/keygen.jl")
include("tfhe/scheme.jl")
export TFHEparams_bin, TFHEparams_block, CCSparams, KMSparams
export CGGI, LMSS, CCS, KMS
export setup, party_keygen, lwe_encrypt, lwe_decrypt, lwe_ith_encrypt, CRS

include("tfhe/bootstrapping.jl")
export bootstrapping!

include("tfhe/gate.jl")
export NAND, AND, OR, XOR, XNOR, NOR, NOT!

include("tfhe/params.jl")
export CGGIparam, Blockparam
export CCS2party, CCS4party, CCS8party, CCS16party
export KMS2party, KMS4party, KMS8party, KMS16party, KMS32party

end