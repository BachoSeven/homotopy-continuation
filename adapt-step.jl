module AdaptStep
  using LinearAlgebra

  export adapt_step

  # Adaptive step size
  function adapt_step(x, x_old, step, m)
    Δ = LinearAlgebra.norm(x - x_old)
    if Δ > 0.1
      step = 0.5 * step
      m = 0
    else
      m+=1
      if (m == 5)
        step = 2 * step
        m = 0
      end
    end

    return (m, step)
  end
end
