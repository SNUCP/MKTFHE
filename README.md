# MKTFHE

<p align="center">
	<img src="logo.png" width="600px"> 
</p>

This implementation is a proof-of-concept for new multi-key TFHE scheme (https://eprint.iacr.org/2022/1460).

For the sake of performance, Float64 is used instead of multi precision floating number when multiplying polynomials with 64 bit coefficients, but the option can be manually turned on when defining the parameters.

Before you run the code, please make sure to install the following packages : ChaChaCiphers, MultiFloats.
To install them, you can open the REPL and type the following commands.

<pre>
<code>
]
add ChaChaCiphers
add MultiFloats
</code>
</pre>

To run the test code for CCS, type the following commands in the terminal.

<pre>
<code>
julia ./test/CCS.jl
</code>
</pre>

To run the test code for our MK-TFHE scheme, type the following commands in the terminal.

<pre>
<code>
julia ./test/KMS.jl
</code>
</pre>

To run the test code for the parallelized version of our MK-TFHE scheme, type the following commands in the terminal.

<pre>
<code>
julia --threads=auto ./test/KMS.jl
</code>
</pre>

This code also provides the implementation of CGGI scheme and LMSS(Faster TFHE bootstrapping from Block Binary Distribution) scheme.
