include("../src/MKTFHE.jl")
using .MKTFHE, Printf

function main()
    params = Blockparam

    @printf("KEY GENERATION ...\n")

    lwekey, ringkey, scheme = setup(params)

    getsize(var) = Base.format_bytes(Base.summarysize(var))
    @printf("BRK SIZE : %s, KSK SIZE : %s\n\n", 
            getsize(scheme.btk.brk), getsize(scheme.btk.ksk))

    gates = [(NAND, (x, y) -> x⊼y, "⊼ "), (AND, (x, y) -> x&y, "& "), (OR, (x, y) -> x||y, "|| "), (XOR, (x, y) -> x⊻y, "⊻ "), (XNOR, (x, y) -> !(x⊻y), "XNOR "), (NOR, (x, y) -> x⊽y, "NOR ")]

    for idx = 1 : 5
        m = rand(Bool, 4)
        ctxts = map(i -> lwe_encrypt(m[i], lwekey, params), 1:4)
        res = ctxts[1]
        mres = m[1]
        circuit = "m1 "

        @time begin
            for i = 2 : 4
                randgate = rand(gates)
                res = randgate[1](res, ctxts[i], scheme)
                mres = randgate[2](mres, m[i])
                circuit *= randgate[3] * "m$i "
            end
        end

        bootstrapping!(res, scheme)
        @assert mres == lwe_decrypt(res, lwekey)
        println("Trial $idx : $circuit= $mres")
    end
end 

main()