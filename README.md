# TFHE

Master branch: [[![CircleCI](https://circleci.com/gh/nucypher/TFHE.jl.svg?style=svg)](https://circleci.com/gh/nucypher/TFHE.jl) [![codecov](https://codecov.io/gh/nucypher/TFHE.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/nucypher/TFHE.jl)
](https://github.com/nucypher/TFHE.jl)

This implementation is a proof-of-concept for new multi-key TFHE scheme (https://eprint.iacr.org/2022/1460).

For the sake of performance, Float64 is used instead of Double64 when multiplying polynomials with 64 bit coefficients, but the option can be manually turned on in the file "polynomials.jl".

Before you run the code, please make sure to install the following packages : GenericFFT, Polynomials, DarkIntegers, DoubleFloats.
To install them, you can open the REPL and type the following commands.

<pre>
<code>
]
add GenericFFT
add DarkIntegers
add DoubleFloats
</code>
</pre>

To run the test code for CCS, type the following commands in the terminal.

<pre>
<code>
julia multikey.jl
</code>
</pre>

To run the test code for Our MK-TFHE scheme, type the following commands in the terminal.

<pre>
<code>
julia multikey_new.jl
</code>
</pre>

To run the test code for the fast implementation of Our MK-TFHE scheme, type the following commands in the terminal.

<pre>
<code>
julia multikey_new_fast.jl
</code>
</pre>
