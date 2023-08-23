module EulerNewton
  using LinearAlgebra
  using TypedPolynomials

  export en_step

  # Euler-Newton predictor-corrector
  function en_step(H, x, t, step_size)

    # Predictor step
    vars = variables(H(t))
    # Jacobian of H evaluated at (x,t)
    JH = [jh(vars=>x) for jh in differentiate(H(t), vars)]
    Δx = JH \ -[gg(vars=>x) for gg in H(1)-H(0)] # ∂H/∂t is the same as γG-F=H(1)-H(0) for our choice of homotopy
    xp = x .+ Δx * step_size

    # Corrector step
    for _ in 1:10
      JH = [jh(vars=>xp) for jh in differentiate(H(t+step_size), vars)]
      Δx = JH \ -[h(vars=>xp) for h in H(t+step_size)]
      xp = xp .+ Δx
      if LinearAlgebra.norm(Δx) < 1e-6
        break
      end
    end

    return xp
  end
end
