# API reference

```@meta
CurrentModule = TFHE
```

## Keys

```@docs
SchemeParameters
tfhe_parameters
SecretKey
CloudKey
make_key_pair
tfhe_parameters_80
tfhe_parameters_128
```

## Encryption/decryption

```@docs
encrypt
decrypt
LweSample
```

## Logical gates

```@docs
gate_nand
gate_or
gate_and
gate_xor
gate_xnor
gate_not
gate_constant
gate_nor
gate_andny
gate_andyn
gate_orny
gate_oryn
gate_mux
```

## Multi-key TFHE (experimental)

```@docs
mktfhe_parameters_2party
mktfhe_parameters_4party
mktfhe_parameters_8party
SharedKey
CloudKeyPart
MKCloudKey
mk_encrypt
mk_decrypt
MKLweSample
mk_gate_nand
```
