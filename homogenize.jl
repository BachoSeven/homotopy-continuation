module Homogenize
  using TypedPolynomials

  export homogenize, homogenized_start_system

  function homogenize(F)
    @polyvar h
    return  [sum([h^(maxdegree(p)-maxdegree(t))*t for t in p.terms]) for p in F]
  end

  function homogenized_start_system(F)
    degrees = [maxdegree(p) for p in F]
    @polyvar h
    G = [x_i^d - h^d for (d, x_i) in zip(degrees, variables(F))]
    r = [[exp(2im*pi/d)^k for k=0:d-1] for d in degrees]
    roots = vec([vcat(collect(root), 1) for root in collect(Iterators.product(r...))])
    return (G, roots)
  end
end
