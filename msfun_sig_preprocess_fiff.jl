using FFTW, Statistics

function msfun_sig_preprocess_fiff(raw::Dict, times::AbstractArray, cfg::Dict)
    if !haskey(cfg, "chans")
        error("cfg must contain a list of channels as 'chans'")
    end

    if typeof(cfg["chans"]) == String
        cfg["chans"] = [cfg["chans"]]
    end

    all_chs = raw["ch_names"]
    chan_inds = findall(ch -> ch in cfg["chans"], all_chs)

    cfg["signal"] = Dict(
        "chan" => chan_inds,
        "ind" => chan_inds,
        "names" => [all_chs[i] for i in chan_inds]
    )

    if isempty(chan_inds)
        println("msfun_sig_preprocess_fiff - WARNING: No channels read... Returning empty output.")
        return zeros(0, size(times, 2)), cfg
    end

    println("msfun_sig_preprocess_fiff - Reading data...")
    sig = raw["data"][chan_inds, :]

    println("msfun_sig_preprocess_fiff - Getting the right time samples...")
    sfreq = raw["sfreq"]
    T = round.(Int, times .* sfreq) .- raw["first_samp"]

    if ndims(times) == 1 || size(times, 1) == 1
        sig = sig[:, T]
    else
        K, L = size(times)
        sig_epo = zeros(K, length(chan_inds), L)
        for k in 1:K
            sig_epo[k, :, :] = sig[:, T[k, :]]
        end
        sig = sig_epo
    end

    if get(cfg, "filter", false)
        println("msfun_sig_preprocess_fiff - Filtering data...")
        if ndims(sig) == 2
            win, F = msfun_prepare_cosine_filter(cfg["filt"], size(sig, 2), sfreq)
            Fsig = fft(sig .* win', 2)
            sig = real.(ifft(Fsig .* F', 2))
        else
            win, F = msfun_prepare_cosine_filter(cfg["filt"], size(sig, 3), sfreq)
            for k in 1:size(sig, 1)
                X = sig[k, :, :]
                FX = fft(X .* win', 2)
                sig[k, :, :] = real.(ifft(FX .* F', 2))
            end
        end
    end

    if get(cfg, "blc", false)
        println("msfun_sig_preprocess_fiff - Applying baseline correction...")
        if ndims(sig) == 2
            n = findall(diff(round.(Int, sfreq .* times)) .> 1)
            n = vcat(0, n, size(times, 2))
            for i in 1:(length(n) - 1)
                avg = mean(sig[:, n[i]+1:n[i+1]], dims=2)
                sig[:, n[i]+1:n[i+1]] .-= avg
            end
        else
            for k in 1:size(sig, 1)
                avg = mean(sig[k, :, :], dims=3)
                sig[k, :, :] .-= avg
            end
        end
    end

    println("msfun_sig_preprocess_fiff - Data preprocessed and ready.")
    return sig, cfg
end