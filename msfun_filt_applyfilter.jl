include("msfun_filt_preparecosine.jl")
using FFTW

function msfun_filt_applyfilter(sig::Array{<:Real}, cfg::Dict)
    if !haskey(cfg, "sfreq") || !haskey(cfg, "filt")
        error("cfg must contain 'sfreq' and 'filt'")
    end

    filt_cfg = cfg["filt"]
    if isa(filt_cfg, String)
        default_filters = Dict(
            "delta"     => Dict("win" => "boxcar", "par" => ["high", "low"], "freq" => [1.0, 4.0], "width" => [0.5, 0.5]),
            "theta"     => Dict("win" => "boxcar", "par" => ["high", "low"], "freq" => [4.0, 8.0], "width" => [1.0, 1.0]),
            "alpha"     => Dict("win" => "boxcar", "par" => ["high", "low"], "freq" => [8.0, 12.0], "width" => [1.0, 1.0]),
            "beta"      => Dict("win" => "boxcar", "par" => ["high", "low"], "freq" => [12.0, 30.0], "width" => [2.0, 2.0]),
            "betalow"   => Dict("win" => "boxcar", "par" => ["high", "low"], "freq" => [12.0, 21.0], "width" => [1.0, 1.0]),
            "betahigh"  => Dict("win" => "boxcar", "par" => ["high", "low"], "freq" => [21.0, 30.0], "width" => [1.0, 1.0]),
            "gamma"     => Dict("win" => "boxcar", "par" => ["high", "low"], "freq" => [30.0, 45.0], "width" => [2.0, 2.0]),
            "gammalow"  => Dict("win" => "boxcar", "par" => ["high", "low"], "freq" => [30.0, 37.5], "width" => [1.0, 1.0]),
            "gammahigh" => Dict("win" => "boxcar", "par" => ["high", "low"], "freq" => [37.5, 45.0], "width" => [1.0, 1.0])
        )
        band = lowercase(filt_cfg)
        if band == "none"
            println("msfun_filt_applyfilter - No filter applied... Just copying data.")
            return copy(sig)
        elseif haskey(default_filters, band)
            filt_cfg = default_filters[band]
        else
            error("Unknown filter name: ", band)
        end
    end

    sfreq = cfg["sfreq"]
    sig_filt = copy(sig)

    println("msfun_filt_applyfilter - Filtering data...")
    if ndims(sig) == 2
        n_ch, T = size(sig)
        win, F = msfun_filt_preparecosine(filt_cfg, T, sfreq)
        win = reshape(win, 1, :)
        F = reshape(F, 1, :)
        Fsig = fft(sig .* win, 2)
        sig_filt = real.(ifft(Fsig .* F, 2))

    elseif ndims(sig) == 3
        n_epochs, n_ch, T = size(sig)
        win, F = msfun_filt_preparecosine(filt_cfg, T, sfreq)
        for k in 1:n_epochs
            epoch = sig[k, :, :]
            Fsig = fft(epoch .* reshape(win, 1, :), 2)
            sig_filt[k, :, :] = real.(ifft(Fsig .* reshape(F, 1, :), 2))
        end
    else
        error("Input must be 2D or 3D")
    end

    println("msfun_filt_applyfilter - Filtered data ready.")
    return sig_filt
end
