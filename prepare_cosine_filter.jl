using FFTW
using DSP

function prepare_cosine_filter(opt::Dict, T::Int, Fs::Real)
    # Apply the specified window function
    win_func = getproperty(DSP, Symbol(opt["win"]))
    win = win_func(T)
    win = collect(transpose(win))  # ensure row vector form

    # Normalize frequency and width by Nyquist (Fs/2)
    norm_freqs = [f / (Fs / 2) for f in opt["freq"]]
    norm_widths = [w / (Fs / 2) for w in opt["width"]]

    # Create cosine filter
    F = cos_filt(T, opt["par"], norm_freqs, norm_widths)
    return win, F
end

function cos_filt(quantum::Int, par_list::Vector{String}, f_vect::Vector, Ws_vect::Vector)
    F_h = ones(quantum)

    for k in eachindex(f_vect)
        f = f_vect[k]
        Ws = Ws_vect[k]
        flow = quantum * (f - Ws / 2) / 2
        fhigh = quantum * (f + Ws / 2) / 2

        F_mult = ones(quantum)

        if par_list[k] == "low"
            F_mult[1:floor(Int, flow)+1] .= 1
            F_mult[ceil(Int, fhigh)+1:ceil(Int, (quantum+1)/2)] .= 0
            trans_idx = ceil(Int, flow)+1:floor(Int, fhigh)+1
            F_mult[trans_idx] = (cos.((trans_idx .- flow .- 1) ./ (fhigh - flow) .* (π / 2))).^2

        elseif par_list[k] == "high"
            F_mult[1:floor(Int, flow)+1] .= 0
            F_mult[ceil(Int, fhigh)+1:ceil(Int, (quantum+1)/2)] .= 1
            trans_idx = ceil(Int, flow)+1:floor(Int, fhigh)+1
            F_mult[trans_idx] = (sin.((trans_idx .- flow .- 1) ./ (fhigh - flow) .* (π / 2))).^2

        elseif par_list[k] == "notch"
            F_mult[1:floor(Int, flow)+1] .= 1
            F_mult[ceil(Int, fhigh)+1:ceil(Int, (quantum+1)/2)] .= 1
            trans_idx = ceil(Int, flow)+1:floor(Int, fhigh)+1
            F_mult[trans_idx] = (cos.((trans_idx .- flow .- 1) ./ (fhigh - flow) .* π)).^2

        else
            error("Filter type must be 'low', 'high', or 'notch'")
        end

        # Mirror for symmetry
        half_idx = ceil(Int, (quantum + 1) / 2)
        F_mult[half_idx+1:end] .= reverse(F_mult[2:half_idx])
        F_h .*= F_mult
    end

    return F_h
end
