# TFHE

Master branch: [![CircleCI](https://circleci.com/gh/nucypher/TFHE.jl.svg?style=svg)](https://circleci.com/gh/nucypher/TFHE.jl) [![codecov](https://codecov.io/gh/nucypher/TFHE.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/nucypher/TFHE.jl)

A Julia port of:

* TFHE (https://github.com/tfhe/tfhe), based on I. Chillotti, N. Gama, M. Georgieva, and M. Izabachène, "Faster Fully Homomorphic Encryption: Bootstrapping in Less Than 0.1 Seconds";
* MKFTHE (https://github.com/ilachill/MK-TFHE), based on H. Chen, I. Chillotti, and Y. Song, "Multi-Key Homomophic Encryption from TFHE".

This implementation is a proof-of-concept for new multi-key TFHE scheme.
For the sake of performance, Float64 is used instead of Double64 when multiplying polynomials with 64 bit coefficients, but the option can be manually turned on in the file "polynomials.jl".

To run the test code for CCS, type julia multikey.jl to the terminal.
To run the test code for Our MK-TFHE scheme, type julia multikey_new.jl to the terminal.
To run the fast implementation of Our MK-TFHE scheme, type julia multikey_new_fast.jl to the terminal.
