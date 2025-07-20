function msfun_filt_concatenate(sig::AbstractArray, L::Int; mode::String = "")
    if ndims(sig) == 2 && isempty(mode)
        mode = "epochlength"
    elseif ndims(sig) == 3 && isempty(mode)
        mode = "epochnum"
    end

    if !(mode in ["epochnum", "epochlength"])
        error("mode must be 'epochnum' or 'epochlength'")
    end

    if ndims(sig) == 2
        N, T = size(sig)
        L = min(L, T)

        if mode == "epochlength"
            epochlength = L
            epochnum = fld(T, L)
        elseif mode == "epochnum"
            epochnum = L
            epochlength = fld(T, L)
        end

        total_len = min(epochnum * epochlength, T)
        sig_trimmed = sig[:, 1:total_len]
        sigbis = reshape(sig_trimmed, N, epochlength, epochnum)
        sigbis = permutedims(sigbis, (3, 1, 2))  # (epochnum, N, epochlength)

    elseif ndims(sig) == 3
        K, N, T = size(sig)
        L = min(L, K)
        sig_trimmed = sig[1:L, :, :]
        sigbis = permutedims(sig_trimmed, (2, 3, 1))  # (N, T, L)
        sigbis = reshape(sigbis, N, L * T)
    else
        error("sig must be 2D or 3D")
    end

    return sigbis
end
