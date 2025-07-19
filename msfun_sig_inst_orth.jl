function msfun_sig_inst_orth(X::AbstractMatrix{ComplexF64}, Y::AbstractMatrix{ComplexF64})
    if size(X, 2) != size(Y, 2)
        error("X and Y must have the same number of time samples (columns)")
    end

    norm_ratio = sum(Y.^2, dims=1) ./ sum(abs.(Y).^2, dims=1)  # shape (1, T)
    proj = conj.(X) .* norm_ratio
    Z = 0.5 .* (X .- proj)

    return Z
end