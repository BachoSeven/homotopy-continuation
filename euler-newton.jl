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
    # ∂H/∂t is the same as γG-F=H(1)-H(0) for our choice of homotopy
    Δx = JH \ -[gg(vars=>x) for gg in H(1)-H(0)]
    xh = x + Δx * step_size

    # Corrector step
    JHh=differentiate(H(t+step_size), vars)
    for _ in 1:10
      JH = [jh(vars=>xh) for jh in JHh]
      Δx = JH \ -[h(vars=>xh) for h in H(t+step_size)]
      xh = xh + Δx
      if LinearAlgebra.norm([h(vars=>xh) for h in H(t+step_size)]) < 1e-8
        break
      end
    end

    return xh
  end
end
