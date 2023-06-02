# TFHE

This implementation is a proof-of-concept for new multi-key TFHE scheme (https://eprint.iacr.org/2022/1460).

For the sake of performance, Float64 is used instead of multi precision floating number when multiplying polynomials with 64 bit coefficients, but the option can be manually turned on when defining the parameters.

Before you run the code, please make sure to install the following packages : Random, Primes, MultiFloats.
To install them, you can open the REPL and type the following commands.

<pre>
<code>
]
add Random
add Primes
add MultiFloats
</code>
</pre>

To run the test code for CCS, type the following commands in the terminal.

<pre>
<code>
julia CCS.jl
</code>
</pre>

To run the test code for Our MK-TFHE scheme, type the following commands in the terminal.

<pre>
<code>
julia KMS.jl
</code>
</pre>

To run the test code for the parallelized version of Our MK-TFHE scheme, type the following commands in the terminal.

<pre>
<code>
julia --threads=auto KMS.jl
</code>
</pre>
