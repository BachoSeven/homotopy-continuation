module EulerNewton
  using LinearAlgebra
  using TypedPolynomials

  export en_step

  # Euler-Newton predictor-corrector
  function en_step(H, x, t, step_size)

    # Predictor step
    vars = variables(H(1))
    # Jacobian of H evaluated at (x,t)
    JH = [jh(vars=>x) for jh in differentiate(H(t), vars)]
    # ∂H/∂t = γG-F = H(1)-H(0) for our homotopy; it doesn't depend on t
    δH_δt = [dh(vars=>x) for dh in H(1)-H(0)]
    Δx = JH \ -δH_δt
    xh = x + Δx * step_size

    # Corrector step
    JHh=differentiate(H(t-step_size), vars)
    for _ in 1:5
      JH = [jh(vars=>xh) for jh in JHh]
      Δx = JH \ -[h(vars=>xh) for h in H(t-step_size)]
      xh = xh + Δx
    end

    return xh
  end
end
