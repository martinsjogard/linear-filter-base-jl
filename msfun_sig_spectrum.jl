using FFTW, Statistics

function msfun_sig_spectrum(sig::AbstractArray, cfg::Dict)
    if !haskey(cfg, "sfreq")
        error("cfg must include key 'sfreq'")
    end

    sfreq = cfg["sfreq"]
    if sfreq <= 0
        error("Sampling frequency must be positive")
    end

    cfg_type = lowercase(get(cfg, "type", "power"))
    if !(cfg_type in ["power", "fourier"])
        error("cfg['type'] must be 'power' or 'fourier'")
    end

    average = get(cfg, "average", false)
    T = size(sig, ndims(sig))
    freq = collect(0:T-1) .* sfreq / T

    Ssig = fft(sig, dims=ndims(sig))
    half = Int(floor(T / 2))
    freq = freq[1:half]

    if ndims(sig) == 2
        Ssig = Ssig[:, 1:half]
    elseif ndims(sig) == 3
        Ssig = Ssig[:, :, 1:half]
    else
        error("Signal must be 2D or 3D")
    end

    if cfg_type == "power"
        Ssig .= abs2.(Ssig)
    end

    if average && ndims(sig) == 3
        if cfg_type == "fourier"
            @warn "msfun_sig_spectrum - Averaging Fourier coefficients is unusual..."
        end
        Ssig = mean(Ssig, dims=1)
        Ssig = dropdims(Ssig, dims=1)
    end

    if get(cfg, "return_band_par", false)
        P = cfg_type == "fourier" ? abs2.(Ssig) : Ssig
        P_sum = sum(P, dims=ndims(P))
        P_norm = P ./ P_sum
        P_cum = cumsum(P_norm, dims=ndims(P))

        band_par = Dict()
        f_idx(x, val) = argmin(abs.(x .- val))
        band_par["fcenter"] = f_idx(P_cum, 0.5)
        band_par["nucenter"] = freq[band_par["fcenter"]]
        band_par["fmin"] = f_idx(P_cum, 1e-2)
        band_par["numin"] = freq[band_par["fmin"]]
        band_par["fmax"] = f_idx(P_cum, 1 - 1e-2)
        band_par["numax"] = freq[band_par["fmax"]]
        return Ssig, freq, band_par
    end

    return Ssig, freq
end