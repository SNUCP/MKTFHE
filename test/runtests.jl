using Base.Iterators: product

using Jute
using Random
using TFHE


gate_tests = [
    ("NAND", gate_nand, 2, !&),
    ("OR", gate_or, 2, |),
    ("AND", gate_and, 2, &),
    ("XOR", gate_xor, 2, xor),
    ("XNOR", gate_xnor, 2, (x, y) -> xor(x, ~y)),
    ("NOT", gate_not, 1, ~),
    ("NOR", gate_nor, 2, !|),
    ("ANDNY", gate_andny, 2, (x, y) -> (~x) & y),
    ("ANDYN", gate_andyn, 2, (x, y) -> x & (~y)),
    ("ORNY", gate_orny, 2, (x, y) -> (~x) | y),
    ("ORYN", gate_oryn, 2, (x, y) -> x | (~y)),
    ("MUX", gate_mux, 3, (x, y, z) -> x ? y : z),
]

gate_test_ids = [gate_test[1] for gate_test in gate_tests]


@testcase "gate" for gate_test in (gate_tests => gate_test_ids)
    rng = MersenneTwister(123)
    secret_key, cloud_key = make_key_pair(rng)

    _, gate, nargs, reference = gate_test

    for bits in product([(false, true) for i in 1:nargs]...)
        ebits = [encrypt(rng, secret_key, b) for b in bits]
        eres = gate(cloud_key, ebits...)
        res = decrypt(secret_key, eres)
        ref_res = reference(bits...)
        @test res == ref_res
    end

end


@testcase "single party, custom parameters" begin
    rng = MersenneTwister(123)
    params = tfhe_parameters_128()
    secret_key, cloud_key = make_key_pair(rng, params)

    _, gate, nargs, reference = gate_tests[1]

    for bits in product([(false, true) for i in 1:nargs]...)
        ebits = [encrypt(rng, secret_key, b) for b in bits]
        eres = gate(cloud_key, ebits...)
        res = decrypt(secret_key, eres)
        ref_res = reference(bits...)
        @test res == ref_res
    end
end


@testcase "multikey NAND" begin

    parties = 2

    params = mktfhe_parameters_2party

    rng = MersenneTwister()

    # Processed on clients' machines
    secret_keys = [SecretKey(rng, params) for i in 1:parties]

    # Created by the server
    shared_key = SharedKey(rng, params)

    # Processed on clients' machines
    ck_parts = [CloudKeyPart(rng, secret_key, shared_key) for secret_key in secret_keys]

    # Processed on the server.
    # `ck_parts` only contain information `public_keys`, `secret_keys` remain secret.
    cloud_key = MKCloudKey(ck_parts)

    for trial = 1:10

        mess1 = rand(Bool)
        mess2 = rand(Bool)
        out = !(mess1 && mess2)

        enc_mess1 = mk_encrypt(rng, secret_keys, mess1)
        enc_mess2 = mk_encrypt(rng, secret_keys, mess2)

        dec_mess1 = mk_decrypt(secret_keys, enc_mess1)
        dec_mess2 = mk_decrypt(secret_keys, enc_mess2)
        @test mess1 == dec_mess1
        @test mess2 == dec_mess2

        enc_out = mk_gate_nand(cloud_key, enc_mess1, enc_mess2)

        dec_out = mk_decrypt(secret_keys, enc_out)
        @test out == dec_out
    end
end


exit(runtests())
