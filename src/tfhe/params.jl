CGGIparam = TFHEparams_bin{UInt32, Float64x2, Float64}(
    630, 2^17,
    8, 2,
    1 << 10, 1, 2^7,
    3, 9
)

Blockparam = TFHEparams_block{UInt32, Float64x2, Float64}(
    329, 3, 2^17,
    8, 2,
    1 << 10, 1, 2^7,
    3, 9
)

CCS2party = CCSparams{UInt32, Float64x2, Float64}(
    560, 130996.5025,
    8, 2,
    1 << 10, 15.9773,
    3, 9,
    2
)

CCS4party = CCSparams{UInt32, Float64x2, Float64}(
    560, 130996.5025,
    8, 2,
    1 << 10, 15.9773,
    4, 8,
    4
)

CCS8party = CCSparams{UInt32, Float64x2, Float64}(
    560, 130996.5025,
    8, 2,
    1 << 10, 15.9773,
    5, 6,
    8
)

CCS16party = CCSparams{UInt32, Float64x2, Float64}(
    560, 130996.5025,
    8, 2,
    1 << 10, 15.9773,
    12, 2,
    16
)

KMS2party = KMSparams{UInt32, UInt64, Float64x2, Float64}(
    560, 130996.5025,
    8, 2,
    1 << 11, 85.4084,
    3, 13, 2, 7, 3, 10,
    2
)

KMS4party = KMSparams{UInt32, UInt64, Float64x2, Float64}(
    560, 130996.5025,
    8, 2,
    1 << 11, 85.4084,
    5, 8, 2, 8, 7, 6,
    4
)

KMS8party = KMSparams{UInt32, UInt64, Float64x2, Float64}(
    560, 130996.5025,
    8, 2,
    1 << 11, 85.4084,
    4, 9, 3, 6, 8, 4,
    8
)

KMS16party = KMSparams{UInt32, UInt64, Float64x2, Float64}(
    560, 130996.5025,
    8, 2,
    1 << 11, 85.4084,
    5, 8, 3, 6, 9, 4,
    16
)

KMS32party = KMSparams{UInt32, UInt64, Float64x2, Float64}(
    560, 130996.5025,
    8, 2,
    1 << 11, 85.4084,
    5, 8, 3, 7, 16, 2,
    32
)