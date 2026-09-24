"""
    LevenbergMarquardt

A trim-safe stand-in for the part of LsqFit.jl that `PositionVelocityTime.user_position` uses:
[`curve_fit`](@ref) with an in-place model and Jacobian,
with the same call signature and the `param` / `resid` fields of its result.

LsqFit keeps the model and Jacobian in abstractly typed fields of NLSolversBase's
`OnceDifferentiable`, so every call through them is a dynamic dispatch, and
`curve_fit` formats caught `MethodError`s and can print a trace — all of which rules
the solve out of a `juliac --trim=safe` build. Here the functions are type
parameters and there is no trace, bounds, timer or error formatting.

The algorithm is LsqFit's `levenberg_marquardt` (after Kanzow, Yamashita & Fukushima
2004) step for step: the same damped normal equations with the `DtD` diagonal
floored at `1e-6`, the same geodesic-acceleration correction when `avv!` is given,
the same step-quality ratio deciding acceptance and trust-region updates, and the
same gradient and relative step-size stopping rules, with the Jacobian evaluated
lazily just as LsqFit's cache does — so it lands on the same iterates.

Kept as a self-contained module, with LsqFit's interface, so that once LsqFit can be
trimmed it can be dropped in favour of `using LsqFit: curve_fit`.
"""
module LevenbergMarquardt

using LinearAlgebra: LAPACK, mul!, norm

"""
    LMResult

The part of `LsqFit.LsqFitResult` that is used: the fitted parameters `param`, the
residual `resid = model(xdata, param) − ydata` at them, and whether the iteration
`converged` (rather than running out of iterations).
"""
struct LMResult
    param::Vector{Float64}
    resid::Vector{Float64}
    converged::Bool
end

"""
    curve_fit(model!, jacobian!, xdata, ydata, p0; inplace = true, avv! = nothing,
              lambda = 10.0, x_tol = 1e-8, g_tol = 1e-12, maxIter = 1000,
              min_step_quality = 1e-3, good_step_quality = 0.75) -> LMResult

Fit `model!(out, xdata, p)` to `ydata` by least squares, starting from `p0`, with
`jacobian!(J, xdata, p)` its Jacobian in `p` — LsqFit's in-place `curve_fit`
(`inplace = true` is the only form provided). `avv!(dir_deriv, p, v)` enables
geodesic acceleration. Like LsqFit, throws an `ArgumentError` if `ydata` holds
non-finite values.
"""
function curve_fit(
    model!::F,
    jacobian!::G,
    xdata,
    ydata::AbstractVector{Float64},
    p0::AbstractVector{Float64};
    inplace::Bool = true,
    avv!::A = nothing,
    lambda::Real = 10.0,
    x_tol::Real = 1e-8,
    g_tol::Real = 1e-12,
    maxIter::Integer = 1000,
    min_step_quality::Real = 1e-3,
    good_step_quality::Real = 0.75,
    lambda_increase::Real = 10.0,
    lambda_decrease::Real = 0.1,
) where {F,G,A}
    inplace || throw(ArgumentError("only the in-place form (`inplace = true`) is provided"))
    all(isfinite, ydata) ||
        throw(ArgumentError("Data contains `Inf` or `NaN` values and a fit cannot be performed"))
    (0 <= min_step_quality < 1) || throw(ArgumentError(" 0 <= min_step_quality < 1 must hold."))
    (0 < good_step_quality <= 1) || throw(ArgumentError(" 0 < good_step_quality <= 1 must hold."))
    (min_step_quality < good_step_quality) ||
        throw(ArgumentError("min_step_quality < good_step_quality must hold."))

    MAX_LAMBDA = 1e16 # minimum trust region radius
    MIN_LAMBDA = 1e-16 # maximum trust region radius
    MIN_DIAGONAL = 1e-6 # lower bound on the diagonal regularising the step

    m = length(ydata)
    n = length(p0)
    x = collect(Float64, p0)
    λ = Float64(lambda)

    # Residual `model − data`, as LsqFit's `curve_fit` forms it.
    residual!(out, p) = (model!(out, xdata, p); out .-= ydata; out)

    f = residual!(Vector{Float64}(undef, m), x)
    J = Matrix{Float64}(undef, m, n)
    jacobian!(J, xdata, x)
    # LsqFit's Jacobian is cached against the point it was last evaluated at and only
    # refreshed at the top of an iteration whose `x` has moved (an accepted step).
    J_is_current = true
    residual = sum(abs2, f)

    trial_f = similar(f)
    trial_x = similar(x)
    JJ = Matrix{Float64}(undef, n, n)
    DtD = Vector{Float64}(undef, n)
    rhs = Vector{Float64}(undef, n)
    v = Vector{Float64}(undef, n)
    a = Vector{Float64}(undef, n)
    delta_x = Vector{Float64}(undef, n)
    Jdelta = Vector{Float64}(undef, m)
    gradient = Vector{Float64}(undef, n)
    dir_deriv = Vector{Float64}(undef, m)

    converged = false
    iteration = 0
    while !converged && iteration < maxIter
        J_is_current || jacobian!(J, xdata, x)
        J_is_current = true

        # Solve (JᵀJ + λ·diag(DtD)) δ = −Jᵀf, with DtD = diag(JᵀJ) floored at
        # MIN_DIAGONAL to prevent "parameter evaporation".
        for i in 1:n
            DtD[i] = max(sum(abs2, view(J, :, i)), MIN_DIAGONAL)
        end
        mul!(JJ, transpose(J), J)
        for i in 1:n
            JJ[i, i] += λ * DtD[i]
        end
        mul!(rhs, transpose(J), f)
        rhs .*= -1
        v .= JJ \ rhs

        if isnothing(avv!)
            delta_x .= v
        else
            # Geodesic acceleration: a = −½ (JᵀJ + λ·diag(DtD))⁻¹ Jᵀ Avv(v).
            avv!(dir_deriv, x, v)
            mul!(a, transpose(J), dir_deriv)
            a .*= -1
            LAPACK.potrf!('U', JJ)
            LAPACK.potrs!('U', JJ, a)
            a .*= 0.5
            delta_x .= v .+ a
        end

        # The residual the linearised model predicts for the step, against the one the
        # step actually achieves.
        mul!(Jdelta, J, delta_x)
        Jdelta .+= f
        predicted_residual = sum(abs2, Jdelta)
        trial_x .= x .+ delta_x
        residual!(trial_f, trial_x)
        trial_residual = sum(abs2, trial_f)

        rho = (trial_residual - residual) / (predicted_residual - residual)
        if trial_residual < residual && rho > min_step_quality
            x .= trial_x
            f .= trial_f
            residual = trial_residual
            J_is_current = false
            if rho > good_step_quality
                λ = max(lambda_decrease * λ, MIN_LAMBDA)
            end
        else
            λ = min(lambda_increase * λ, MAX_LAMBDA)
        end
        iteration += 1

        # As in LsqFit, the gradient pairs the (not yet refreshed) Jacobian of this
        # iteration with the residual after its step.
        mul!(gradient, transpose(J), f)
        converged =
            norm(gradient, Inf) < g_tol || norm(delta_x) < x_tol * (x_tol + norm(x))
    end
    return LMResult(x, f, converged)
end

end
