function msfun_filt_preparecosine(opt::Dict, T::Int, Fs::Float64)
    # Generate time-domain window
    win = ones(Float64, T)
    if opt["win"] == "hann"
        win .= 0.5 .- 0.5 * cos.(2π .* collect(0:T-1) ./ (T - 1))
    elseif opt["win"] == "boxcar"
        win .= 1.0
    else
        error("Unsupported window type: $(opt["win"])")
    end

    # Normalize frequencies by Nyquist (Fs / 2)
    norm_freqs = opt["freq"] ./ Fs .* 2
    norm_widths = opt["width"] ./ Fs .* 2

    # Compute cosine filter
    F = cos_filt(T, opt["par"], norm_freqs, norm_widths)
    return win, F
end

function cos_filt(quantum::Int, par_list::Vector{String}, f_vect::Vector{Float64}, Ws_vect::Vector{Float64})
    F_h = ones(Float64, quantum)

    for (k, (f, Ws)) in enumerate(zip(f_vect, Ws_vect))
        flow = quantum * (f - Ws / 2) / 2
        fhigh = quantum * (f + Ws / 2) / 2

        F_mult = ones(Float64, quantum)

        if par_list[k] == "low"
            F_mult[1:floor(Int, flow)+1] .= 1
            F_mult[ceil(Int, fhigh)+1:ceil(Int, (quantum+1)/2)] .= 0
            trans_idx = ceil(Int, flow)+1:floor(Int, fhigh)+1
            F_mult[trans_idx] .= (cos.((trans_idx .- flow .- 1) ./ (fhigh - flow) .* (π / 2))).^2

        elseif par_list[k] == "high"
            F_mult[1:floor(Int, flow)+1] .= 0
            F_mult[ceil(Int, fhigh)+1:ceil(Int, (quantum+1)/2)] .= 1
            trans_idx = ceil(Int, flow)+1:floor(Int, fhigh)+1
            F_mult[trans_idx] .= (sin.((trans_idx .- flow .- 1) ./ (fhigh - flow) .* (π / 2))).^2

        elseif par_list[k] == "notch"
            F_mult[1:floor(Int, flow)+1] .= 1
            F_mult[ceil(Int, fhigh)+1:ceil(Int, (quantum+1)/2)] .= 1
            trans_idx = ceil(Int, flow)+1:floor(Int, fhigh)+1
            F_mult[trans_idx] .= (cos.((trans_idx .- flow .- 1) ./ (fhigh - flow) .* π)).^2

        else
            error("The filter must be 'low', 'high', or 'notch'")
        end

        half_idx = ceil(Int, (quantum + 1) / 2)
        F_mult[half_idx+1:end] .= reverse(F_mult[2:half_idx])
        F_h .*= F_mult
    end

    return F_h
end
