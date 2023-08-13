"""
bootstrapping!(ctxt, scheme) bootstraps a binary BFV ciphertext.
"""
function bootstrapping!(ctxt::LWE{T}, scheme::TFHEscheme{R, S}) where {T, R, S}
    N, logN = scheme.N, trailing_zeros(scheme.N)

    # Scale up the TLWE ciphertext to a LWE ciphertext with modulus 2N.
    tildea = divbits.(ctxt.a, bits(T) - logN - 1)
    tildeb = divbits(ctxt.b, bits(T) - logN - 1)

    oneovereight = R(1) << (bits(R) - 3)
    b = zeronativepoly(N, R)
    if tildeb ≤ N
        @inbounds @simd for i = 1 : N
            b.coeffs[i] = i ≤ tildeb ? oneovereight : -oneovereight
        end
    else
        tildeb -= R(N)
        @inbounds @simd for i = 1 : N
            b.coeffs[i] = i ≤ tildeb ? -oneovereight : oneovereight
        end
    end
    acc = RLWE(b, [zeronativepoly(N, R) for _ = 1 : scheme.k])

    blindrotate!(tildea, acc, scheme)
    keyswitch!(ctxt, acc, scheme)
end

"""
Blind Rotation algorithm based on CGGI16.
"""
function blindrotate!(tildea::Vector{<:Unsigned}, acc::RLWE{R}, scheme::CGGI{R, S}) where {R, S}
    N = scheme.ffter.N

    param = scheme.gswpar

    # pre-allocate the vector to store the result of RLWE decomposition.
    bvec = [buffnativepoly(N, R) for _ = 1 : param.l]
    avec = [buffnativepoly(N, R) for _ = 1 : param.l, _ = 1 : param.k]
    tbvec = [bufftransnativepoly(N, S) for _ = 1 : param.l]
    tavec = [bufftransnativepoly(N, S) for _ = 1 : param.l, _ = 1 : param.k]

    # pre-allocate the temporary accumulator.
    acc2 = deepcopy(acc)
    tacc = fft(acc, scheme.ffter)

    for idx = eachindex(tildea)
        if tildea[idx] > 0 
            # decompose the accumulator RLWE ciphertext.
            decompto!(bvec, acc.b, param)
            decompto!(avec, acc.a, param)

            # transform the decomposed polynomial into evaluation form.
            @inbounds @simd for i = eachindex(bvec)
                fftto!(tbvec[i], bvec[i], scheme.ffter)
            end
            @inbounds @simd for i = eachindex(avec)
                fftto!(tavec[i], avec[i], scheme.ffter)
            end
            
            # perform external product.
            initialise!(tacc)
            @inbounds for i = 1 : param.l
                muladdto!(tacc, tbvec[i], scheme.btk.brk[idx].basketb.stack[i])
            end
            @inbounds for i = 1 : param.k, j = 1 : param.l
                muladdto!(tacc, tavec[j, i], scheme.btk.brk[idx].basketa[i].stack[j])
            end
            
            # acc = acc + (X^ai - 1) acc2
            mul!(scheme.monomial[tildea[idx]], tacc)
            ifftto!(acc2, tacc, scheme.ffter)
            add!(acc, acc2)
        end
    end
end

"""
Key-switch algorithm based on CGGI16.
"""
function keyswitch!(res::LWE{T}, acc::RLWE{T}, scheme::CGGI{T, S}) where {T, S}
    n, N, k, l = scheme.n, scheme.N, scheme.k, scheme.kskpar.l

    # initialise the resulting ciphertext and set b.
    initialise!(res)
    res.b = acc.b.coeffs[1]

    # iteratively perform key-switching, with on-the-fly ciphertext extraction.
    ajdec = Vector{T}(undef, scheme.kskpar.l)
    @inbounds @simd for i = 1 : k
        unbalanceddecompto!(ajdec, acc.a[i].coeffs[1], scheme.kskpar)
        @inbounds @simd for j = 1 : l
            if ajdec[j] > 0
                add!(res, scheme.btk.ksk[ajdec[j], 1, i].stack[j])
            end
        end

        @inbounds @simd for idx = 2 : N
            unbalanceddecompto!(ajdec, T(0) - acc.a[i].coeffs[N - idx + 2], scheme.kskpar)
            @inbounds @simd for j = 1 : l
                if ajdec[j] > 0
                    add!(res, scheme.btk.ksk[ajdec[j], idx, i].stack[j])
                end
            end
        end
    end

    res
end

"""
Blind Rotation algorithm based on LMSS23 (Fast TFHE Bootstrapping from Block Binary Distribution).
"""
function blindrotate!(tildea::Vector{<:Unsigned}, acc::RLWE{R}, scheme::LMSS{R, S}) where {R, S}
    N = scheme.ffter.N

    # pre-allocate the vector to store decomposed polynomials.
    param = scheme.gswpar
    bvec = [buffnativepoly(N, R) for _ = 1 : param.l]
    avec = [buffnativepoly(N, R) for _ = 1 : param.l, _ = 1 : param.k]
    tbvec = [bufftransnativepoly(N, S) for _ = 1 : param.l]
    tavec = [bufftransnativepoly(N, S) for _ = 1 : param.l, _ = 1 : param.k]

    # pre-allocate the temporary accumulators.
    acc2 = deepcopy(acc)
    tacc = fft(acc, scheme.ffter)
    tacc2 = deepcopy(tacc)

    for idx1 = 0 : scheme.d-1
        # decompose the accumulator.
        decompto!(bvec, acc.b, param)
        decompto!(avec, acc.a, param)

        # transform the decomposed polynomials into evaluation form.
        @inbounds @simd for i = eachindex(bvec)
            fftto!(tbvec[i], bvec[i], scheme.ffter)
        end
        @inbounds @simd for i = eachindex(avec)
            fftto!(tavec[i], avec[i], scheme.ffter)
        end

        initialise!(tacc2)
        for idx2 = 1 : scheme.ℓ
            idx = idx1 * scheme.ℓ + idx2
            if tildea[idx] > 0
                initialise!(tacc)

                # perform external product, and store the result in tacc.
                @inbounds for i = 1 : param.l
                    muladdto!(tacc, tbvec[i], scheme.btk.brk[idx].basketb.stack[i])
                end
                @inbounds for i = 1 : param.k, j = 1 : param.l
                    muladdto!(tacc, tavec[j, i], scheme.btk.brk[idx].basketa[i].stack[j])
                end
                
                # add to tacc2, with multiplying (X^aᵢ - 1)
                muladdto!(tacc2, scheme.monomial[tildea[idx]], tacc)
            end
        end
        
        # transform tacc2 back to coefficient form and add to accumulator.
        ifftto!(acc2, tacc2, scheme.ffter)
        add!(acc, acc2)
    end
end

"""
Key-switch algorithm based on LMSS23.
"""
function keyswitch!(res::LWE{T}, acc::RLWE{T}, scheme::LMSS{T, S}) where {T, S}
    n, N, k, l = scheme.n, scheme.N, scheme.k, scheme.kskpar.l

    initialise!(res)
    res.b = acc.b.coeffs[1]
    current = 1

    ajdec = Vector{T}(undef, scheme.kskpar.l)
    for i = 1 : k
        if current + N ≤ n 
            res.a[current] = acc.a[i].coeffs[1]
            @inbounds @simd for j = 1 : N-1
                res.a[current+j] = T(0) - acc.a[i].coeffs[N-j+1]
            end
            current += N

            continue
        elseif current ≤ n 
            res.a[current] = acc.a[i].coeffs[1]
            @inbounds @simd for j = 1 : n - current
                res.a[current+j] = T(0) - acc.a[i].coeffs[N-j+1]
            end

            @inbounds @simd for idx = n-current+2 : N
                decompto!(ajdec, T(0) - acc.a[i].coeffs[N - idx + 2], scheme.kskpar)
                @inbounds @simd for j = 1 : l
                    if signed(ajdec[j]) > 0
                        add!(res, scheme.btk.ksk[ajdec[j], idx, i].stack[j])
                    elseif ajdec[j] != 0
                        sub!(res, scheme.btk.ksk[-signed(ajdec[j]), idx, i].stack[j])
                    end
                end
            end

            current = n+1
        else
            decompto!(ajdec, acc.a[i].coeffs[1], scheme.kskpar)
            @inbounds @simd for j = 1 : l
                if signed(ajdec[j]) > 0
                    add!(res, scheme.btk.ksk[ajdec[j], 1, i].stack[j])
                elseif ajdec[j] != 0
                    sub!(res, scheme.btk.ksk[-signed(ajdec[j]), 1, i].stack[j])
                end
            end

            @inbounds @simd for idx = 2 : N
                decompto!(ajdec, T(0) - acc.a[i].coeffs[N - idx + 2], scheme.kskpar)
                @inbounds @simd for j = 1 : l
                    if signed(ajdec[j]) > 0
                        add!(res, scheme.btk.ksk[ajdec[j], idx, i].stack[j])
                    elseif ajdec[j] != 0
                        sub!(res, scheme.btk.ksk[-signed(ajdec[j]), idx, i].stack[j])
                    end
                end
            end
        end
    end

    res
end

"""
Blind Rotation algorithm based on CCS19.
"""
function blindrotate!(tildeavec::Vector{<:Unsigned}, acc::RLWE{T}, scheme::CCS{T, S}) where {T, S}
    tildea = reshape(tildeavec, scheme.n, scheme.k)

    k, n, N = scheme.k, scheme.n, scheme.N
   
    param = scheme.unipar

    bvec = [buffnativepoly(N, T) for _ = 1 : param.l]
    avec = [buffnativepoly(N, T) for _ = 1 : param.l, _ = 1 : k]
    tbvec = [bufftransnativepoly(N, S) for _ = 1 : param.l]
    tavec = [bufftransnativepoly(N, S) for _ = 1 : param.l, _ = 1 : k]

    v0 = buffnativepoly(N, T)
    tv0 = bufftransnativepoly(N, S)
    v0vec = [buffnativepoly(N, T) for _ = 1 : param.l]
    vvec = [buffnativepoly(N, T) for _ = 1 : param.l, _ = 1 : k]

    v = [buffnativepoly(N, T) for _ = 1 : k]
    tv = [bufftransnativepoly(N, S) for _ = 1 : k]
    tv0vec = [bufftransnativepoly(N, S) for _ = 1 : param.l]
    tvvec = [bufftransnativepoly(N, S) for _ = 1 : param.l, _ = 1 : k]

    tacc = TransRLWE(bufftransnativepoly(N, S), [bufftransnativepoly(N, S) for _ = 1 : k])
    acc2 = RLWE(buffnativepoly(N, T), [buffnativepoly(N, T) for _ = 1 : k])

    @inbounds for idx = 1 : k
        @inbounds for i = 1 : n
            if tildea[i, idx] > 0
                # hybrid product
                # decompose b and a of acc.
                decompto!(bvec, acc.b, param)
                decomptoith!(avec, acc.a, idx, param)

                # Transform bvec and avec to FFT form.
                @inbounds @simd for j = 1 : param.l
                    fftto!(tbvec[j], bvec[j], scheme.ffter)
                end
                @inbounds @simd for j = 1 : idx
                    @inbounds @simd for k = 1 : param.l 
                        fftto!(tavec[k, j], avec[k, j], scheme.ffter)
                    end
                end

                # Compute u.
                initialise!(tacc)
                @inbounds for j = 1 : param.l
                    muladdto!(tacc.b, tbvec[j], scheme.btk[idx].brk[i].d[j])
                end
                @inbounds for j = 1 : idx, k = 1 : param.l
                    muladdto!(tacc.a[j], tavec[k, j], scheme.btk[idx].brk[i].d[k])
                end

                # Compute v.
                @. initialise!(tv)
                initialise!(tv0)
                @inbounds for j = 1 : param.l
                    mulsubto!(tv0, tbvec[j], scheme.a[j])
                end
                @inbounds for j = 1 : idx, k = 1 : param.l
                    muladdto!(tv[j], tavec[k, j], scheme.btk[j].b[k])
                end

                # Inverse transform tv to v.
                ifftto!(v0, tv0, scheme.ffter)
                @inbounds for j = 1 : idx
                    ifftto!(v[j], tv[j], scheme.ffter)
                end

                # Decompose v and transfomr.
                decompto!(v0vec, v0, param)
                decomptoith!(vvec, v, idx, param)
                @inbounds for j = 1 : param.l
                    fftto!(tv0vec[j], v0vec[j], scheme.ffter)
                end 
                @inbounds for j = 1 : idx, k = 1 : param.l
                    fftto!(tvvec[k, j], vvec[k, j], scheme.ffter)
                end

                # Compute w.
                @inbounds for j = 1 : param.l
                    muladdto!(tacc.b, tv0vec[j], scheme.btk[idx].brk[i].f.stack[j].b)
                    muladdto!(tacc.a[idx], tv0vec[j], scheme.btk[idx].brk[i].f.stack[j].a[1])
                end
                @inbounds for j = 1 : idx, k = 1 : param.l
                    muladdto!(tacc.b, tvvec[k, j], scheme.btk[idx].brk[i].f.stack[k].b)
                    muladdto!(tacc.a[idx], tvvec[k, j], scheme.btk[idx].brk[i].f.stack[k].a[1])
                end

                mul!(scheme.monomial[tildea[i, idx]], tacc)
                ifftto!(acc2, tacc, scheme.ffter)
                add!(acc, acc2)
            end
        end
    end
end

"""
Key-switch algorithm based on CCS19.
"""
function keyswitch!(res::LWE{T}, acc::RLWE{T}, scheme::CCS{T, S}) where {T, S}
    n, N, k, l = scheme.n, scheme.N, scheme.k, scheme.kskpar.l

    initialise!(res)
    res.b = acc.b.coeffs[1]

    # defining partial LWE ciphertext reduces the allocations.
    partctxt = [LWE(zero(T), zeros(T, n)) for _ = 1 : k]

    # key-switching is parallelizable.
    @threads for i = 1 : k
        ajdec = Vector{T}(undef, l)
        unbalanceddecompto!(ajdec, acc.a[i].coeffs[1], scheme.kskpar)
        @inbounds @simd for idx = eachindex(ajdec)
            if ajdec[idx] > 0
                add!(partctxt[i], scheme.btk[i].ksk[ajdec[idx], 1].stack[idx])
            end
        end

        @inbounds @simd for j = 2 : N
            unbalanceddecompto!(ajdec, T(0) - acc.a[i].coeffs[N - j + 2], scheme.kskpar)
            @inbounds for idx = eachindex(ajdec)
                if ajdec[idx] > 0
                    add!(partctxt[i], scheme.btk[i].ksk[ajdec[idx], j].stack[idx])
                end
            end
        end

        res.b += partctxt[i].b
        @. res.a[(i-1)*n+1:i*n] += partctxt[i].a
    end
end

"""
Blind Rotation algorithm based on our scheme.
"""
function blindrotate!(tildeavec::Vector{T}, acc::RLWE{R}, scheme::KMSScheme{T, R, S}) where {T, R, S}
    tildea = reshape(tildeavec, scheme.n, scheme.k)
    
    k = scheme.k
    levkey = Vector{TransRLEV{S}}(undef, k)

    # phase 1 : party-wise single key blind rotation.
    @sync for i = 1 : k
        @spawn levkey[i] = phase_1(i, tildea[:, i], scheme.btk[i], scheme.monomial, scheme.ffter)
    end

    # phase 2 : merge them into a single MK ciphertext.
    phase_2!(levkey, acc, scheme)

    acc
end

"""
Phase 1 from Blind Rotation of our scheme.
"""
function phase_1(party::Int64, tildea::Vector{T}, btk::BootKey_KMS{T, R, S}, monomial::Vector{TransNativePoly{S}}, ffter::FFTransformer{S}) where {T, R, S}
    N = ffter.N

    bvec = [buffnativepoly(N, R) for _ = 1 : btk.gswpar.l]
    avec = [buffnativepoly(N, R) for _ = 1 : btk.gswpar.l]
    tbvec = [bufftransnativepoly(N, S) for _ = 1 : btk.gswpar.l]
    tavec = [bufftransnativepoly(N, S) for _ = 1 : btk.gswpar.l]

    # For party 1, we don't need to perform blind rotation on RLEV ciphertext.
    # This implementation may lead to a bigger noise variance, rather than performing blind rotation with test vector.
    # We implement phase 1 this way to keep the code simpler.
    iter = party == 1 ? 1 : btk.levpar.l

    stack = Vector{RLWE{R}}(undef, iter)
    @inbounds @simd for i = 1 : iter
        stack[i] = RLWE(zeronativepoly(N, R), zeronativepoly(N, R))
        stack[i].b.coeffs[1] = btk.levpar.gvec[i]
    end
    acc = RLEV(stack)
    acc2 = RLEV(deepcopy(stack))
    tacc = fft(acc, ffter)

    # Basically the same code to blind rotation of CGGI, but for each row of the RLEV ciphertext.
    @inbounds for idx = eachindex(tildea)
        if tildea[idx] > 0 
            @inbounds @simd for i = 1 : iter
                decompto!(bvec, acc.stack[i].b, btk.gswpar)
                decompto!(avec, acc.stack[i].a[1], btk.gswpar)

                initialise!(tacc.stack[i])

                @inbounds @simd for j = eachindex(bvec)
                    fftto!(tbvec[j], bvec[j], ffter)
                end
                @inbounds @simd for j = eachindex(avec)
                    fftto!(tavec[j], avec[j], ffter)
                end

                @inbounds for j = eachindex(bvec)
                    muladdto!(tacc.stack[i], tbvec[j], btk.brk[idx].basketb.stack[j])
                end
                @inbounds for j = eachindex(avec)
                    muladdto!(tacc.stack[i], tavec[j], btk.brk[idx].basketa[1].stack[j])
                end
            end

            mul!(monomial[tildea[idx]], tacc)
            ifftto!(acc2, tacc, ffter)
            add!(acc, acc2)
        end
    end

    fftto!(tacc, acc, ffter)
    tacc
end

"""
Phase 2 from Blind Rotation of our scheme.
"""
function phase_2!(levkey::Vector{TransRLEV{S}}, acc::RLWE{R}, scheme::KMSScheme{T, R, S}) where {T, R, S}
    k, N = scheme.k, scheme.N
    
    levpar, unipar = scheme.btk[1].levpar, scheme.btk[1].unipar
    maxl = max(levpar.l, unipar.l)

    bvec = [buffnativepoly(N, R) for _ = 1 : maxl]
    avec = [buffnativepoly(N, R) for _ = 1 : maxl, _ = 1 : k]
    tbvec = [bufftransnativepoly(N, S) for _ = 1 : maxl]
    tavec = [bufftransnativepoly(N, S) for _ = 1 : maxl, _ = 1 : k]

    v = buffnativepoly(N, R)
    tv = bufftransnativepoly(N, S)
    vvec = [buffnativepoly(N, R) for _ = 1 : unipar.l]
    tvvec = [bufftransnativepoly(N, S) for _ = 1 : unipar.l]

    tx = TransRLWE(bufftransnativepoly(N, S), [bufftransnativepoly(N, S) for _ = 1 : k])
    y = RLWE(buffnativepoly(N, R), [buffnativepoly(N, R) for _ = 1 : k])
    ty = TransRLWE(bufftransnativepoly(N, S), [bufftransnativepoly(N, S) for _ = 1 : k])

    @inbounds for idx = 1 : k
        # decompose b and ai of acc
        decompto!(bvec, acc.b, levpar)
        decomptoith!(avec, acc.a, idx-1, levpar)

        # FFT the decomposed bvec and avec
        @inbounds for i = 1 : levpar.l
            fftto!(tbvec[i], bvec[i], scheme.ffter)
        end
        @inbounds for i = 1 : idx-1, j = 1 : levpar.l
            fftto!(tavec[j, i], avec[j, i], scheme.ffter)
        end

        iter = idx == 1 ? 1 : levpar.l

        # LEV multiplication between bvec and levkey[i]
        initialise!(tx)
        @inbounds for i = 1 : iter
            muladdto!(tx.b, tbvec[i], levkey[idx].stack[i].b)
        end
        @inbounds for i = 1 : idx-1, j = 1 : iter
            muladdto!(tx.a[i], tavec[j, i], levkey[idx].stack[j].b)
        end

        # LEV multiplication between avec and levkey[i]
        initialise!(ty)
        @inbounds for i = 1 : iter
            muladdto!(ty.b, tbvec[i], levkey[idx].stack[i].a[1])
        end
        @inbounds for i = 1 : idx-1, j = 1 : iter
            muladdto!(ty.a[i], tavec[j, i], levkey[idx].stack[j].a[1])
        end

        ifftto!(y.b, ty.b, scheme.ffter)
        @inbounds for i = 1 : idx-1
            ifftto!(y.a[i], ty.a[i], scheme.ffter)
        end

        # hybrid product
        # decompose b and a of yi
        decompto!(bvec, y.b, unipar)
        decomptoith!(avec, y.a, idx-1, unipar)

        # Transform bvec and avec in FFT form
        @inbounds for i = 1 : unipar.l
            fftto!(tbvec[i], bvec[i], scheme.ffter)
        end
        @inbounds for i = 1 : idx-1, j = 1 : unipar.l
            fftto!(tavec[j, i], avec[j, i], scheme.ffter)
        end

        # Compute u. Reuse ty for memory save.
        initialise!(ty)
        @inbounds for i = 1 : unipar.l
            muladdto!(ty.b, tbvec[i], scheme.btk[idx].rlk.d[i])
        end
        @inbounds for i = 1 : idx-1, j = 1 : unipar.l
            muladdto!(ty.a[i], tavec[j, i], scheme.btk[idx].rlk.d[j])
        end

        # Compute v.
        initialise!(tv)
        @inbounds for i = 1 : unipar.l
            mulsubto!(tv, tbvec[i], scheme.a[i])
        end
        @inbounds for i = 1 : idx - 1, j = 1 : unipar.l
            muladdto!(tv, tavec[j, i], scheme.btk[i].b[j])
        end

        # Inverse transform tv to v.
        ifftto!(v, tv, scheme.ffter)

        # Decompose v and trasform.
        decompto!(vvec, v, unipar)
        @inbounds for i = 1 : unipar.l
            fftto!(tvvec[i], vvec[i], scheme.ffter)
        end

        # Compute w.
        @inbounds for i = 1 : unipar.l
            muladdto!(ty.b, tvvec[i], scheme.btk[idx].rlk.f.stack[i].b)
            muladdto!(ty.a[idx], tvvec[i], scheme.btk[idx].rlk.f.stack[i].a[1])
        end

        # add back to x.
        add!(tx, ty)

        # Inverse transform to acc.
        ifftto!(acc, tx, scheme.ffter)
    end 
end

"""
Key-switch algorithm from our scheme.
Basically the same code to CCS19, but with modulus change from R to T.
"""
function keyswitch!(res::LWE{T}, acc::RLWE{R}, scheme::KMS{T, R, S}) where {T, R, S}
    n, N, k, l = scheme.n, scheme.N, scheme.k, scheme.kskpar.l

    bitdiff = bits(R) - bits(T)
    initialise!(res)
    res.b = T(acc.b.coeffs[1] >> bitdiff)

    partctxt = [LWE(zero(T), zeros(T, n)) for _ = 1 : k]

    @threads for i = 1 : k
        ajdec = Vector{T}(undef, l)
        unbalanceddecompto!(ajdec, T(acc.a[i].coeffs[1] >> bitdiff), scheme.kskpar)
        @inbounds @simd for idx = eachindex(ajdec)
            if ajdec[idx] > 0
                add!(partctxt[i], scheme.btk[i].ksk[ajdec[idx], 1].stack[idx])
            end
        end

        @inbounds @simd for j = 2 : N
            unbalanceddecompto!(ajdec, T(0) - T(acc.a[i].coeffs[N - j + 2] >> bitdiff), scheme.kskpar)
            @inbounds @simd for idx = eachindex(ajdec)
                if ajdec[idx] > 0
                    add!(partctxt[i], scheme.btk[i].ksk[ajdec[idx], j].stack[idx])
                end
            end
        end

        res.b += partctxt[i].b
        @. res.a[(i-1)*n+1:i*n] += partctxt[i].a
    end
end

"""
Phase 1 from Blind Rotation of KMS with block binary keys.
"""
function phase_1(party::Int64, tildea::Vector{T}, btk::BootKey_KMS_block{T, R, S}, monomial::Vector{TransNativePoly{S}}, ffter::FFTransformer{S}) where {T, R, S}
    N = ffter.N

    bvec = [buffnativepoly(N, R) for _ = 1 : btk.gswpar.l]
    avec = [buffnativepoly(N, R) for _ = 1 : btk.gswpar.l]
    tbvec = [bufftransnativepoly(N, S) for _ = 1 : btk.gswpar.l]
    tavec = [bufftransnativepoly(N, S) for _ = 1 : btk.gswpar.l]

    # For party 1, we don't need to perform blind rotation on RLEV ciphertext.
    # This implementation may lead to a bigger noise variance, rather than performing blind rotation with test vector.
    # We implement phase 1 this way to keep the code simpler.
    iter = party == 1 ? 1 : btk.levpar.l

    stack = Vector{RLWE{R}}(undef, iter)
    @inbounds @simd for i = 1 : iter
        stack[i] = RLWE(zeronativepoly(N, R), zeronativepoly(N, R))
        stack[i].b.coeffs[1] = btk.levpar.gvec[i]
    end
    acc = RLEV(stack)
    acc2 = RLEV(deepcopy(stack))
    tacc = fft(acc, ffter)
    tacc2 = deepcopy(tacc)

    # Basically the same code to blind rotation of LMSS, but for each row of the RLEV ciphertext.
    @inbounds for idx1 = 0 : btk.d-1
        @inbounds @simd for i = 1 : iter
            decompto!(bvec, acc.stack[i].b, btk.gswpar)
            decompto!(avec, acc.stack[i].a[1], btk.gswpar)

            @inbounds @simd for j = eachindex(bvec)
                fftto!(tbvec[j], bvec[j], ffter)
            end
            @inbounds @simd for j = eachindex(avec)
                fftto!(tavec[j], avec[j], ffter)
            end

            initialise!(tacc2.stack[i])
            for idx2 = 1 : btk.ℓ
                idx = idx1 * btk.ℓ + idx2
                if tildea[idx] > 0     
                    initialise!(tacc.stack[i])

                    @inbounds for j = eachindex(bvec)
                        muladdto!(tacc.stack[i], tbvec[j], btk.brk[idx].basketb.stack[j])
                    end
                    @inbounds for j = eachindex(avec)
                        muladdto!(tacc.stack[i], tavec[j], btk.brk[idx].basketa[1].stack[j])
                    end

                    muladdto!(tacc2.stack[i], monomial[tildea[idx]], tacc.stack[i])
                end
            end
        end

        ifftto!(acc2, tacc2, ffter)
        add!(acc, acc2)
    end

    fftto!(tacc, acc, ffter)
    tacc
end

"""
Key-switch algorithm for KMS with block binary keys.
"""
function keyswitch!(res::LWE{T}, acc::RLWE{R}, scheme::KMS_block{T, R, S}) where {T, R, S}
    n, N, k, l = scheme.n, scheme.N, scheme.k, scheme.kskpar.l

    bitdiff = bits(R) - bits(T)
    initialise!(res)
    res.b = T(acc.b.coeffs[1] >> bitdiff)

    partctxt = [LWE(zero(T), zeros(T, n)) for _ = 1 : k]

    @threads for i = 1 : k
        ajdec = Vector{T}(undef, l)

        partctxt[i].a[1] = T(acc.a[i].coeffs[1] >> bitdiff)
        @inbounds @simd for j = 2 : n
            partctxt[i].a[j] = T(0) - T(acc.a[i].coeffs[N-j+2] >> bitdiff)
        end

        @inbounds @simd for j = n+1 : N
            decompto!(ajdec, T(0) - T(acc.a[i].coeffs[N - j + 2] >> bitdiff), scheme.kskpar)
            @inbounds @simd for idx = eachindex(ajdec)
                if signed(ajdec[idx]) > 0
                    add!(partctxt[i], scheme.btk[i].ksk[ajdec[idx], j].stack[idx])
                elseif ajdec[idx] != 0
                    sub!(partctxt[i], scheme.btk[i].ksk[-signed(ajdec[idx]), j].stack[idx])
                end
            end
        end

        res.b += partctxt[i].b
        @. res.a[(i-1)*n+1:i*n] += partctxt[i].a
    end
end