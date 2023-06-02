function NAND(ctxt1::LWE{T}, ctxt2::LWE{T}, scheme::TFHEscheme) where T
    b = T(1) << (bits(T) - 3) - ctxt1.b - ctxt2.b
    a = zeros(T, length(ctxt1.a))
    @. a = - ctxt1.a - ctxt2.a
    res = LWE(b, a)
    bootstrapping!(res, scheme)
    res
end

function AND(ctxt1::LWE{T}, ctxt2::LWE{T}, scheme::TFHEscheme) where T
    b = T(7) << (bits(T) - 3) + ctxt1.b + ctxt2.b
    a = zeros(T, length(ctxt1.a))
    @. a = ctxt1.a + ctxt2.a
    res = LWE(b, a)
    bootstrapping!(res, scheme)
    res
end

function OR(ctxt1::LWE{T}, ctxt2::LWE{T}, scheme::TFHEscheme) where T
    b = T(1) << (bits(T) - 3) + ctxt1.b + ctxt2.b
    a = zeros(T, length(ctxt1.a))
    @. a = ctxt1.a + ctxt2.a
    res = LWE(b, a)
    bootstrapping!(res, scheme)
    res
end

function XOR(ctxt1::LWE{T}, ctxt2::LWE{T}, scheme::TFHEscheme) where T
    b = T(1) << (bits(T) - 2) + T(2) * (ctxt1.b + ctxt2.b)
    a = zeros(T, length(ctxt1.a))
    @. a = T(2) * (ctxt1.a + ctxt2.a)
    res = LWE(b, a)
    bootstrapping!(res, scheme)
    res
end

function XNOR(ctxt1::LWE{T}, ctxt2::LWE{T}, scheme::TFHEscheme) where T
    b = T(3) << (bits(T) - 2) - T(2) * (ctxt1.b + ctxt2.b)
    a = zeros(T, length(ctxt1.a))
    @. a = -T(2) * (ctxt1.a + ctxt2.a)
    res = LWE(b, a)
    bootstrapping!(res, scheme)
    res
end

function NOR(ctxt1::LWE{T}, ctxt2::LWE{T}, scheme::TFHEscheme) where T
    b = T(7) << (bits(T) - 3) - ctxt1.b - ctxt2.b
    a = zeros(T, length(ctxt1.a))
    @. a = - ctxt1.a - ctxt2.a
    res = LWE(b, a)
    bootstrapping!(res, scheme)
    res
end

function NOT!(ctxt::LWE{T}) where T
    ctxt.b *= -T(1)
    @. ctxt.a *= -T(1)
end