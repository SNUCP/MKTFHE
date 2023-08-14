CGGIparam = TFHEparams_bin{UInt32, Float64x2, Float64}(
    630, 1 << 17,
    8, 2,
    1 << 10, 1, 1 << 7,
    3, 9
)

Blockparam = TFHEparams_block{UInt32, Float64x2, Float64}(
    229, 3, 1 << 17,
    8, 2,
    1 << 10, 1, 1 << 7,
    3, 9
)

CCS2party = CCSparams{UInt32, Float64x2, Float64}(
    560, 1 << 17,
    8, 2,
    1 << 10, 1 << 4,
    3, 8,
    2
)

CCS4party = CCSparams{UInt32, Float64x2, Float64}(
    560, 1 << 17,
    8, 2,
    1 << 10, 1 << 4,
    4, 8,
    4
)

CCS8party = CCSparams{UInt32, Float64x2, Float64}(
    560, 1 << 17,
    8, 2,
    1 << 10, 1 << 4,
    5, 6,
    8
)

CCS16party = CCSparams{UInt32, Float64x2, Float64}(
    560, 1 << 17,
    8, 2,
    1 << 10, 1 << 4,
    12, 2,
    16
)

KMS2party = KMSparams{UInt32, UInt64, Float64x2, Float64}(
    560, 1 << 17,
    8, 2,
    1 << 11, 85.4084,
    3, 12, 2, 7, 3, 10,
    2
)

KMS4party = KMSparams{UInt32, UInt64, Float64x2, Float64}(
    560, 1 << 17,
    8, 2,
    1 << 11, 85.4084,
    5, 8, 2, 8, 7, 6,
    4
)

KMS8party = KMSparams{UInt32, UInt64, Float64x2, Float64}(
    560, 1 << 17,
    8, 2,
    1 << 11, 85.4084,
    4, 9, 3, 6, 8, 4,
    8
)

KMS16party = KMSparams{UInt32, UInt64, Float64x2, Float64}(
    560, 1 << 17,
    8, 2,
    1 << 11, 85.4084,
    5, 8, 3, 6, 9, 4,
    16
)

KMS32party = KMSparams{UInt32, UInt64, Float64x2, Float64}(
    560, 1 << 17,
    8, 2,
    1 << 11, 85.4084,
    6, 7, 3, 7, 16, 2,
    32
)

KMS2partyblock = KMSparams_block{UInt32, UInt64, Float64x2, Float64}(
    203, 3, 1 << 17,
    8, 2,
    1 << 11, 85.4084,
    3, 12, 2, 7, 3, 10,
    2
)

KMS4partyblock = KMSparams_block{UInt32, UInt64, Float64x2, Float64}(
    203, 3, 1 << 17,
    8, 2,
    1 << 11, 85.4084,
    5, 8, 2, 8, 7, 6,
    4
)

KMS8partyblock = KMSparams_block{UInt32, UInt64, Float64x2, Float64}(
    203, 3, 1 << 17,
    8, 2,
    1 << 11, 85.4084,
    4, 9, 3, 6, 8, 4,
    8
)

KMS16partyblock = KMSparams_block{UInt32, UInt64, Float64x2, Float64}(
    203, 3, 1 << 17,
    8, 2,
    1 << 11, 85.4084,
    5, 8, 3, 6, 9, 4,
    16
)

KMS32partyblock = KMSparams_block{UInt32, UInt64, Float64x2, Float64}(
    203, 3, 1 << 17,
    8, 2,
    1 << 11, 85.4084,
    6, 7, 3, 7, 16, 2,
    32
)