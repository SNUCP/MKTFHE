abstract type BootKey end

struct BootKey_bin{T, S} <: BootKey where {T<:Unsigned, S<:AbstractFloat}
    brk::Vector{TransRGSW{S}}
    ksk::Array{LEV{T}, 3}

    function BootKey_bin(lwekey::LWEkey{T}, ringkey::RLWEkey{T}, kskpar::LEVparams_digit{T}, α::Float64,
                         gswpar::GSWparams_digit{T}, β::Float64, ffter_keygen::FFTransformer{R}, ffter::FFTransformer{S}) where {T, R, S}
        n, k, N, D = lwekey.n, ringkey.k, ringkey.N, 1 << kskpar.logB
        brk = Vector{TransRGSW{S}}(undef, n)
        ksk = Array{LEV{T}, 3}(undef, D-1, N, k)
    
        @threads for i = 1 : n
            brk[i] = fft(rgsw_encrypt(lwekey.key[i], ringkey, β, ffter_keygen, gswpar), ffter)
        end
    
        @inbounds @simd for idx = 1 : k
            @threads for i = 1 : N
                @inbounds @simd for j = 1 : D-1
                    ksk[j, i, idx] = lev_encrypt(T(ringkey.key[idx].coeffs[i] * j), lwekey, α, kskpar)
                end
            end
        end
    
        new{T, S}(brk, ksk)
    end
end

struct BootKey_block{T, S} <: BootKey where {T<:Unsigned, S<:AbstractFloat}
    brk::Vector{TransRGSW{S}}
    ksk::Array{LEV{T}, 3}

    function BootKey_block(lwekey::LWEkey{T}, ringkey::RLWEkey{T}, kskpar::LEVparams_digit{T}, α::Float64,
                           gswpar::GSWparams_digit{T}, β::Float64, ffter_keygen::FFTransformer{R}, ffter::FFTransformer{S}) where {T, R, S}
        n, k, N, D = lwekey.n, ringkey.k, ringkey.N, 1 << kskpar.logB
        brk = Vector{TransRGSW{S}}(undef, n)
        ksk = Array{LEV{T}, 3}(undef, D÷2, N, k)
    
        @threads for i = 1 : n
            brk[i] = fft(rgsw_encrypt(lwekey.key[i], ringkey, β, ffter_keygen, gswpar), ffter)
        end
    
        @inbounds @simd for idx = 1 : k
            @threads for i = 1 : N
                if (idx-1)*N+i > n
                    @inbounds @simd for j = 1 : D÷2
                        ksk[j, i, idx] = lev_encrypt(T(ringkey.key[idx].coeffs[i] * j), lwekey, α, kskpar)
                    end
                end
            end
        end
    
        new{T, S}(brk, ksk)
    end
end

struct BootKey_CCS{T, S} <: BootKey where {T<:Unsigned, S<:AbstractFloat}
    b::Vector{TransNativePoly{S}}
    brk::Vector{TransUniEnc{S}}
    ksk::Array{LEV{T}, 2}

    function BootKey_CCS(lwekey::LWEkey{T}, ringkey::RLWEkey{T}, kskpar::LEVparams_digit{T}, α::Float64, 
                         a::CRS{T}, unipar::Uniparams_digit{T}, β::Float64, ffter_keygen::FFTransformer{R}, ffter::FFTransformer{S}) where {T, R, S}
        n, N, D = lwekey.n, ringkey.N, 1 << kskpar.logB

        ta = fft.(a, Ref(ffter_keygen))
        b = fft.(gen_b(ta, ringkey, β, ffter_keygen, unipar), Ref(ffter))
        brk = Vector{TransUniEnc{S}}(undef, n)
        ksk = Array{LEV{T}, 2}(undef, D-1, N)
    
        @threads for i = 1 : n
            brk[i] = fft(unienc_encrypt(ta, lwekey.key[i], ringkey, β, ffter_keygen, unipar), ffter)
        end
    
        @threads for i = 1 : N
            @inbounds @simd for j = 1 : D-1
                ksk[j, i] = lev_encrypt(T(ringkey.key[1].coeffs[i] * j), lwekey, α, kskpar)
            end
        end
    
        new{T, S}(b, brk, ksk)
    end
end

struct BootKey_KMS{T, R, S} <: BootKey where {T<:Unsigned, R<:Unsigned, S<:AbstractFloat}
    b::Vector{TransNativePoly{S}}
    brk::Vector{TransRGSW{S}}
    rlk::TransUniEnc{S}
    ksk::Array{LEV{T}, 2}
    gswpar::GSWparams_digit{R}
    levpar::LEVparams_digit{R}
    unipar::Uniparams_digit{R}

    function BootKey_KMS(lwekey::LWEkey{T}, gswkey::RLWEkey{R}, unikey::RLWEkey{R}, kskpar::LEVparams_digit, α::Float64, 
                         levpar::LEVparams_digit{R}, gswpar::GSWparams_digit{R}, β::Float64, a::CRS, unipar::Uniparams_digit{R}, 
                         ffter_keygen::FFTransformer{U}, ffter::FFTransformer{S}) where {T, R, U, S}
        n, N, D = lwekey.n, gswkey.N, 1 << kskpar.logB

        ta = fft.(a, Ref(ffter_keygen))
        b = fft.(gen_b(ta, unikey, β, ffter_keygen, unipar), Ref(ffter))

        brk = Vector{TransRGSW{S}}(undef, n)
        rlk = fft(unienc_encrypt(ta, gswkey.key[1], unikey, β, ffter_keygen, unipar), ffter)
        ksk = Array{LEV{T}, 2}(undef, D-1, N)
    
        @threads for i = 1 : n
            brk[i] = fft(rgsw_encrypt(R(lwekey.key[i]), gswkey, β, ffter_keygen, gswpar), ffter)
        end

        @threads for i = 1 : N
            @inbounds @simd for j = 1 : D-1
                ksk[j, i] = lev_encrypt(T(unikey.key[1].coeffs[i] * j), lwekey, α, kskpar)
            end
        end
    
        new{T, R, S}(b, brk, rlk, ksk, gswpar, levpar, unipar)
    end
end