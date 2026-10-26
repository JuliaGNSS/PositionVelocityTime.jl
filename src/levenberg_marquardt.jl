"""
    LevenbergMarquardt

A trim-safe stand-in for the part of LsqFit.jl that `PositionVelocityTime.estimate_position`
uses: LsqFit's in-place `curve_fit`, as [`curve_fit!`](@ref) on a reusable
[`LMWorkspace`](@ref), with LsqFit's arguments after the workspace and the `param` /
`resid` fields of its result.

LsqFit calls the model and Jacobian through abstractly typed fields of NLSolversBase's
`OnceDifferentiable` and formats caught `MethodError`s, which a `juliac --trim=safe` build
rejects. Here the functions are type parameters, with no trace, bounds, timer or error
formatting.

The algorithm is LsqFit's `levenberg_marquardt` (after Kanzow, Yamashita & Fukushima
2004) step for step, down to its lazily refreshed Jacobian, so it lands on the same
iterates to rounding (the damped normal equations are solved by Cholesky rather than LU).
On a reusable workspace a fit allocates nothing.

Once LsqFit can be trimmed, this module can be dropped for `using LsqFit: curve_fit`:
`curve_fit!(LMWorkspace(), args...)` is that call.
"""
module LevenbergMarquardt

using LinearAlgebra: LAPACK, mul!, norm

"""
    LMResult

The part of `LsqFit.LsqFitResult` that is used: the fitted parameters `param`, the
residual `resid = model(xdata, param) − ydata` at them, and whether the iteration
`converged` before running out of iterations.
"""
struct LMResult
    param::Vector{Float64}
    resid::Vector{Float64}
    converged::Bool
end

"""
    LMWorkspace()

The scratch storage of [`curve_fit!`](@ref): every vector and matrix of the iteration,
grown to the largest problem seen and reused. The matrices, which cannot be resized in
place, are used through views of their leading `m × n` block.
"""
mutable struct LMWorkspace
    x::Vector{Float64}
    f::Vector{Float64}
    trial_f::Vector{Float64}
    trial_x::Vector{Float64}
    DtD::Vector{Float64}
    rhs::Vector{Float64}
    v::Vector{Float64}
    a::Vector{Float64}
    delta_x::Vector{Float64}
    Jdelta::Vector{Float64}
    gradient::Vector{Float64}
    dir_deriv::Vector{Float64}
    J::Matrix{Float64}
    JJ::Matrix{Float64}
end

LMWorkspace() = LMWorkspace(
    Float64[], Float64[], Float64[], Float64[], Float64[], Float64[],
    Float64[], Float64[], Float64[], Float64[], Float64[], Float64[],
    Matrix{Float64}(undef, 0, 0), Matrix{Float64}(undef, 0, 0))

"""
    grown(matrix, rows, cols) -> Matrix{Float64}

`matrix` if it holds at least `rows × cols` elements in that shape, else a new matrix
that does. The caller stores the result back and works on
`view(matrix, 1:rows, 1:cols)`.
"""
grown(matrix::Matrix{Float64}, rows, cols) =
    size(matrix, 1) >= rows && size(matrix, 2) >= cols ? matrix :
    Matrix{Float64}(undef, max(rows, size(matrix, 1)), max(cols, size(matrix, 2)))

"""
    curve_fit!(workspace::LMWorkspace, model!, jacobian!, xdata, ydata, p0; kwargs...)
        -> LMResult

Fit `model!(out, xdata, p)` to `ydata` by least squares from `p0`, with
`jacobian!(J, xdata, p)` its Jacobian (handed a view, not a `Matrix`) — LsqFit's in-place
`curve_fit` on the buffers of `workspace`. Allocates nothing once they have grown.
`avv!(dir_deriv, p, v)` enables geodesic acceleration. Throws an `ArgumentError` for
non-finite `ydata`, as LsqFit does.

The keywords and their defaults are LsqFit's: `inplace = true` (the only form provided),
`avv! = nothing`, `lambda = 10.0`, `x_tol = 1e-8`, `g_tol = 1e-12`, `maxIter = 1000`,
`min_step_quality = 1e-3`, `good_step_quality = 0.75`, `lambda_increase = 10.0` and
`lambda_decrease = 0.1`.

The result's `param` and `resid` **are** the workspace's `x` and `f`, overwritten by the
next fit; `p0` may itself be `workspace.x`.
"""
function curve_fit!(
    workspace::LMWorkspace,
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
    # `copyto!` onto itself is a no-op, which is what makes `p0 === workspace.x` safe.
    x = copyto!(resize!(workspace.x, n), p0)
    λ = Float64(lambda)

    # Residual `model − data`, as LsqFit's `curve_fit` forms it.
    residual!(out, p) = (model!(out, xdata, p); out .-= ydata; out)

    workspace.J = grown(workspace.J, m, n)
    workspace.JJ = grown(workspace.JJ, n, n)
    J = view(workspace.J, 1:m, 1:n)
    JJ = view(workspace.JJ, 1:n, 1:n)

    f = residual!(resize!(workspace.f, m), x)
    jacobian!(J, xdata, x)
    # As LsqFit's cache: refreshed only at the top of an iteration whose `x` has moved.
    J_is_current = true
    residual = sum(abs2, f)

    trial_f = resize!(workspace.trial_f, m)
    trial_x = resize!(workspace.trial_x, n)
    DtD = resize!(workspace.DtD, n)
    rhs = resize!(workspace.rhs, n)
    v = resize!(workspace.v, n)
    a = resize!(workspace.a, n)
    delta_x = resize!(workspace.delta_x, n)
    Jdelta = resize!(workspace.Jdelta, m)
    gradient = resize!(workspace.gradient, n)
    dir_deriv = resize!(workspace.dir_deriv, m)

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
        # The damped normal matrix is positive definite, so one Cholesky factor serves the
        # step and the geodesic correction (LsqFit's `\` allocates an LU per iteration;
        # the two agree to rounding). Should it fail anyway — the damping lost in
        # rounding — the step counts as rejected: damp harder and try again.
        _, info = LAPACK.potrf!('U', JJ)
        if info != 0
            λ = min(lambda_increase * λ, MAX_LAMBDA)
            iteration += 1
            continue
        end
        LAPACK.potrs!('U', JJ, copyto!(v, rhs))

        if isnothing(avv!)
            delta_x .= v
        else
            # Geodesic acceleration: a = −½ (JᵀJ + λ·diag(DtD))⁻¹ Jᵀ Avv(v).
            avv!(dir_deriv, x, v)
            mul!(a, transpose(J), dir_deriv)
            a .*= -1
            LAPACK.potrs!('U', JJ, a)
            a .*= 0.5
            delta_x .= v .+ a
        end

        # The residual the linearised model predicts for the step, and the one it achieves.
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

        # As in LsqFit: this iteration's Jacobian with the residual after its step.
        mul!(gradient, transpose(J), f)
        converged =
            norm(gradient, Inf) < g_tol || norm(delta_x) < x_tol * (x_tol + norm(x))
    end
    return LMResult(x, f, converged)
end

end
