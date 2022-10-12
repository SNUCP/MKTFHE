# Manual

For more detailed examples see `examples/tutorial.jl`, `examples/multikey.jl` and `test/runtests.jl`.


## A simple example

TFHE works on separate bits and has all the common logical gates defined for encrypted bits: AND, OR, NAND and so on.
Each gate is bootstrapped, so the error of any ciphertext is always small enough for it to be decrypted, and bootstrapping never has to be called explicitly.

In this simple example we will create two vectors of encrypted bits and run them through a XOR gate.
First, we create the secret key (used for encryption and decryption) and the cloud key (as follows from the name, used for secure computation by a third party).

```julia
using TFHE
using Random

rng = MersenneTwister(123)
secret_key, cloud_key = make_key_pair(rng)
```

Then we encrypt two short vectors of booleans.

```julia
bits1 = [false, true, false, true]
bits2 = [false, false, true, true]

ciphertext1 = encrypt.(Ref(rng), secret_key, bits1)
ciphertext2 = encrypt.(Ref(rng), secret_key, bits2)
```

XOR gate is applied on encrypted bits using the cloud key.

```julia
cresult = gate_xor.(cloud_key, ciphertext1, ciphertext2)

```

Calculate the reference result and check that the decryption coincides with it.

```julia
reference = xor.(bits1, bits2)
result = decrypt.(secret_key, cresult)
@assert result == reference
```
