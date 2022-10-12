include("./src/TFHE.jl")
using Random, Printf
using .TFHE


function main()
    rng = MersenneTwister()
    test(rng)
end

main()