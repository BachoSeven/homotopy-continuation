module AdaptStep
  using LinearAlgebra
  using TypedPolynomials

  export adapt_step

  # Adaptive step size
  function adapt_step(H, x, t, step, m)
    Δ = norm([h(variables(H(t))=>x) for h in H(t-step)])
    if Δ > 1e-10
      step = 0.5 * step
      m = 0
    else
      m+=1
      if (m == 4)
        step = 2 * step
        m = 0
      end
    end

    return (m, step)
  end
end
