using Statistics

function msfun_sig_downsample(sig::AbstractArray, cfg::Dict)
    if ndims(sig) == 2
        K = 1
        C, T = size(sig)
    elseif ndims(sig) == 3
        K, C, T = size(sig)
    else
        error("sig must be 2D or 3D")
    end

    sfreq = cfg["sfreq"]
    downsfreq = cfg["downsfreq"]
    smooth = get(cfg, "smooth", true)
    overlap = get(cfg, "overlap", 1)

    N = sfreq / downsfreq
    if N % 1 != 0
        println("msfun_sig_downsample - WARNING: sfreq/downsfreq not integer, rounding")
        N = round(Int, N)
        downsfreq = sfreq / N
        println("New downsampling frequency: ", downsfreq)
    end

    N = Int(N)
    if overlap >= N
        overlap = 1
    end

    step = fld(N, overlap)
    nsamp = fld(T - N, step) + 1
    tsamp = round(Int, N / 2) .+ (0:step:(nsamp - 1) * step)

    if ndims(sig) == 3
        sig_flat = reshape(permutedims(sig, (3, 1, 2)), T, K * C)'
    else
        sig_flat = sig
    end

    if !smooth
        sigbis = sig_flat[:, tsamp]
    else
        tbuf = [i + (0:N-1) for i in 0:step:(nsamp - 1) * step]
        sigbuf = hcat([sig_flat[:, idx] for idx in tbuf]...)
        sigbuf = reshape(sigbuf, size(sig_flat, 1), N, nsamp)
        sigbis = dropdims(mean(sigbuf, dims=2), dims=2)
    end

    if ndims(sig) == 3
        sigbis = reshape(sigbis', nsamp, K, C)
        sigbis = permutedims(sigbis, (2, 3, 1))
    end

    return sigbis, tsamp
end