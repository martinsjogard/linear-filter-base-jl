using FFTW

function msfun_sig_slow_modulation(Z::AbstractArray, fcenter::Int)
    if ndims(Z) != 2 && ndims(Z) != 3
        error("Z must be a 2D or 3D array")
    end

    T = size(Z, ndims(Z))
    isanalytic = eltype(Z) <: Complex

    if !isanalytic
        Z = msfun_sig_analytic(Z, ndims(Z))
    end

    phase_wave = exp.(2π * im * (fcenter - 1) * (0:T-1) ./ T)
    shape = size(Z)
    temp_shape = ntuple(i -> i == ndims(Z) ? 1 : shape[i], ndims(Z))
    phase_wave = reshape(phase_wave, temp_shape...)

    Zslow = Z ./ phase_wave

    if !isanalytic
        Zslow = real(Zslow)
    end

    return Zslow
end