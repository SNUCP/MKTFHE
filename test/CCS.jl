using MKTFHE, Printf

function main()
    params = CCS2party
    a = CRS(params)

    @printf("KEY GENERATION ...\n")

    keys = [party_keygen(a, params) for _ = 1 : params.k]
    lwekeys = first.(keys)
    btk = last.(keys)

    getsize(var) = Base.format_bytes(Base.summarysize(var))
    @printf("BRK SIZE : %s, KSK SIZE : %s\n\n", 
            getsize(btk[1].brk), getsize(btk[1].ksk))

    scheme = setup(a, btk, params)

    gates = [(NAND, (x, y) -> x⊼y, "⊼ "), (AND, (x, y) -> x&y, "& "), (OR, (x, y) -> x||y, "|| "), (XOR, (x, y) -> x⊻y, "⊻ "), (XNOR, (x, y) -> !(x⊻y), "XNOR "), (NOR, (x, y) -> x⊽y, "NOR ")]

    for idx = 1 : 5
        m = rand(Bool, params.k)
        ctxts = map(i -> lwe_ith_encrypt(m[i], i, lwekeys[i], params), 1 : params.k)
        res = ctxts[1]
        mres = m[1]
        circuit = "m1 "

        for i = 2 : params.k
            randgate = rand(gates)
            res = randgate[1](res, ctxts[i], scheme)
            mres = randgate[2](mres, m[i])
            circuit *= randgate[3] * "m$i "
        end

        @time bootstrapping!(res, scheme)
        @assert mres == lwe_decrypt(res, lwekeys, params)
        println("Trial $idx : $circuit= $mres")
    end
end

main()