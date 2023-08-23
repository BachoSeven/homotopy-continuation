module Homotopy
  export homotopy

  # Define a straight-line homotopy between the two systems
  function homotopy(F, G)
    γ = cis(2π * rand())
    function H(t)
      return [(1 - t) * f + γ * t * g for (f, g) in zip(F, G)]
    end
    return H
  end
end
