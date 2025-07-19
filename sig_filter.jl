using FFTW
include("prepare_cosine_filter.jl")

function sig_filter(sig::Array, cfg::Dict)
    # Check signal dimensions
    nd = ndims(sig)
    if nd != 2 && nd != 3
        error("sig must be 2D or 3D")
    end

    if !haskey(cfg, "sfreq") || !haskey(cfg, "filt")
        error("cfg must contain 'sfreq' and 'filt'")
    end

    filt_cfg = cfg["filt"]

    # Handle default named filters
    if isa(filt_cfg, String)
        name = lowercase(filt_cfg)
        if name == "none"
            println("sig_filter - WARNING : No filter applied... Just copying data.")
            return copy(sig)
        end
        default_filters = Dict(
            "delta" => Dict("win" => "rect", "par" => ["high", "low"], "freq" => [1.0, 4.0], "width" => [0.5, 0.5]),
            "theta" => Dict("win" => "rect", "par" => ["high", "low"], "freq" => [4.0, 8.0], "width" => [1.0, 1.0]),
            "alpha" => Dict("win" => "rect", "par" => ["high", "low"], "freq" => [8.0, 12.0], "width" => [1.0, 1.0]),
            "beta" => Dict("win" => "rect", "par" => ["high", "low"], "freq" => [12.0, 30.0], "width" => [2.0, 2.0]),
            "betalow" => Dict("win" => "rect", "par" => ["high", "low"], "freq" => [12.0, 21.0], "width" => [1.0, 1.0]),
            "betahigh" => Dict("win" => "rect", "par" => ["high", "low"], "freq" => [21.0, 30.0], "width" => [1.0, 1.0]),
            "gamma" => Dict("win" => "rect", "par" => ["high", "low"], "freq" => [30.0, 45.0], "width" => [2.0, 2.0]),
            "gammalow" => Dict("win" => "rect", "par" => ["high", "low"], "freq" => [30.0, 37.5], "width" => [1.0, 1.0]),
            "gammahigh" => Dict("win" => "rect", "par" => ["high", "low"], "freq" => [37.5, 45.0], "width" => [1.0, 1.0])
        )
        if haskey(default_filters, name)
            cfg["filt"] = default_filters[name]
        else
            error("Unknown filter name: ", name)
        end
    end

    # Copy input signal
    sig_filt = copy(sig)
    sfreq = cfg["sfreq"]
    filt_cfg = cfg["filt"]

    println("sig_filter - Filtering data...")
    if nd == 2
        # Continuous
        win, F = prepare_cosine_filter(filt_cfg, size(sig, 2), sfreq)
        win = reshape(win, 1, :)
        F = reshape(F, 1, :)
        Fsig = fft(sig .* win, 2)
        sig_filt = real.(ifft(Fsig .* F, 2))
    else
        # Epoched
        win, F = prepare_cosine_filter(filt_cfg, size(sig, 3), sfreq)
        for k in 1:size(sig, 1)
            epoch = sig[k, :, :]
            Fsig = fft(epoch .* reshape(win, 1, :), 2)
            sig_filt[k, :, :] = real.(ifft(Fsig .* reshape(F, 1, :), 2))
        end
    end

    return sig_filt
end
