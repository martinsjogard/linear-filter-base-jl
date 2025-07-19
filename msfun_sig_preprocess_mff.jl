using FFTW, Statistics

function msfun_sig_preprocess_mff(times::AbstractArray, sfreq::Real, cfg::Dict)
    if !haskey(cfg, "chans") || !haskey(cfg, "mff_data")
        error("cfg must include 'mff_data' and 'chans'")
    end

    data = cfg["mff_data"]
    if typeof(cfg["chans"]) == String
        cfg["chans"] = [cfg["chans"]]
    end

    chan_indices = findall(ch -> ch in cfg["chans"], cfg["channel_names"])
    if isempty(chan_indices)
        println("No valid channels found, returning empty output.")
        return zeros(0, size(times, 2)), cfg
    end

    sig = data[chan_indices, :]
    cfg["signal"] = Dict(
        "chan" => chan_indices,
        "names" => cfg["chans"]
    )

    Tidx = round.(Int, times .* sfreq)
    Tidx .-= minimum(Tidx)

    if ndims(times) == 1 || size(times, 1) == 1
        sig = sig[:, Tidx]
    else
        K, L = size(times)
        sig_epoch = zeros(K, length(chan_indices), L)
        for k in 1:K
            sig_epoch[k, :, :] = sig[:, Tidx[k, :]]
        end
        sig = sig_epoch
    end

    if get(cfg, "filter", false)
        println("Applying filter...")
        win, F = msfun_prepare_cosine_filter(cfg["filt"], size(sig, ndims(sig)), sfreq)
        if ndims(sig) == 2
            Fsig = fft(sig .* win', 2)
            sig = real.(ifft(Fsig .* F', 2))
        else
            for k in 1:size(sig, 1)
                X = sig[k, :, :]
                FX = fft(X .* win', 2)
                sig[k, :, :] = real.(ifft(FX .* F', 2))
            end
        end
    end

    if get(cfg, "blc", false)
        println("Applying baseline correction...")
        if ndims(sig) == 2
            n = findall(diff(round.(Int, sfreq .* times)) .> 1)
            pushfirst!(n, 0)
            push!(n, size(times, 2))
            for i in 1:(length(n)-1)
                avg = mean(sig[:, n[i]+1:n[i+1]], dims=2)
                sig[:, n[i]+1:n[i+1]] .-= avg
            end
        else
            for k in 1:size(sig, 1)
                avg = mean(sig[k, :, :], dims=3)
                for t in 1:size(sig, 3)
                    sig[k, :, t] .-= avg
                end
            end
        end
    end

    println("msfun_sig_preprocess_mff - Data preprocessed and ready.")
    return sig, cfg
end