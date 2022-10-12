include("./src/TFHE.jl")
using Random, Printf
using .TFHE

function main()
    parties = 2

    params = mktfhe_parameters_2party

    rng = MersenneTwister()

    @printf("KEY GENERATION ...\n")

    function keygen()
        # Processed on clients' machines
        secret_keys = [SecretKey(rng, params) for _ in 1:parties]

        # Created by the server
        shared_key = SharedKey(rng, params)

        # Processed on clients' machines
        ck_parts = [CloudKeyPart(rng, secret_key, shared_key) for secret_key in secret_keys]

        # Processed on the server.
        # `ck_parts` only contain information `public_keys`, `secret_keys` remain secret.
        secret_keys, MKCloudKey(ck_parts, shared_key)
    end

    @time secret_keys, cloud_key = keygen()

    getsize(var) = Base.format_bytes(Base.summarysize(var)/parties)

    @printf("BK SIZE : %s, KSK SIZE : %s\n\n", 
            getsize(cloud_key.bootstrap_key), getsize(cloud_key.keyswitch_key))

    for trial = 1:5

        mess1 = rand(Bool)
        mess2 = rand(Bool)
        out = !(mess1 && mess2)

        enc_mess1 = mk_encrypt(rng, secret_keys, mess1)
        enc_mess2 = mk_encrypt(rng, secret_keys, mess2)

        dec_mess1 = mk_decrypt(secret_keys, enc_mess1)
        dec_mess2 = mk_decrypt(secret_keys, enc_mess2)
        @assert mess1 == dec_mess1
        @assert mess2 == dec_mess2

        @time enc_out = mk_gate_nand(cloud_key, enc_mess1, enc_mess2)

        dec_out = mk_decrypt(secret_keys, enc_out)
        @assert out == dec_out

        println("Trial $trial: $mess1 NAND $mess2 = $out")
    end
end


main()