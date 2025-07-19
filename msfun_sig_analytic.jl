using FFTW

function msfun_sig_analytic(X::AbstractArray, dim::Int = ndims(X))
    if !(eltype(X) <: Real)
        @warn "msfun_sig_analytic - Input is not real. Using real part only."
        X = real(X)
    end

    if dim < 1 || dim > ndims(X)
        error("Invalid dimension")
    end

    if size(X, dim) == 1
        error("Selected dimension contains only one sample")
    end

    X_perm = permutedims(X, (dim, setdiff(1:ndims(X), dim)...))
    n = size(X_perm, 1)

    Xf = fft(X_perm, 1)
    Zf = zeros(ComplexF64, size(Xf))

    if iseven(n)
        Zf[1, :] .= Xf[1, :]
        Zf[2:n÷2, :] .= 2 .* Xf[2:n÷2, :]
        Zf[n÷2+1, :] .= Xf[n÷2+1, :]
    else
        Zf[1, :] .= Xf[1, :]
        Zf[2:(n+1)÷2, :] .= 2 .* Xf[2:(n+1)÷2, :]
    end

    Z = ifft(Zf, 1)
    invperm = invpermute((dim, setdiff(1:ndims(X), dim)...))
    return permutedims(Z, invperm)
end