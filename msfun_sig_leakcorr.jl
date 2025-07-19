function msfun_sig_leakcorr(X::AbstractArray, Y::AbstractArray, cfg::Dict)
    if ndims(X) != 2 || ndims(Y) != 2 || size(X, 2) != size(Y, 2)
        error("X and Y must be 2D arrays with matching time dimension")
    end

    method = lowercase(get(cfg, "method", ""))
    if !(method in ["gcs", "orthinst", "orthstat", "custom"])
        error("cfg['method'] must be one of 'gcs', 'orthinst', 'orthstat', 'custom'")
    end

    if method == "gcs"
        if size(Y, 1) != 1
            error("GCS method requires Y to have shape (1, T)")
        end
        gcs = cfg["gcs"]
        inv = gcs["inv"]
        ind = gcs["ind"]
        beta = inv["invop"] * inv["leadfield"][:, ind]
        beta ./= beta[ind]
        Z = X .- beta * Y[1, :]

    elseif method == "orthstat"
        beta = (real(X) * real(Y)') * pinv(real(Y) * real(Y)')
        Z = X .- beta * Y

    elseif method == "orthinst"
        if eltype(X) <: Real && eltype(Y) <: Real
            error("X and Y must be complex for 'orthinst'")
        end
        ratio = sum(Y .^ 2, dims=1) ./ sum(abs.(Y) .^ 2, dims=1)
        Z = 0.5 .* (X .- conj.(X) .* ratio)

    elseif method == "custom"
        beta = cfg["beta"]
        if size(beta) != (size(X, 1), size(Y, 1))
            error("cfg['beta'] must be of shape (N, M)")
        end
        Z = X .- beta * Y
    end

    return Z
end