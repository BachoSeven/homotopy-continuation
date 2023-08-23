module StartSystem
  using TypedPolynomials

  export start_system

  # Define start system based on total degree
  function start_system(F)
    degrees = [maxdegree(p) for p in F]
    G = [x_i^d - 1 for (d, x_i) in zip(degrees, variables(F))]
    r = [[exp(2im*pi/d)^k for k=0:d-1] for d in degrees]
    roots = vec([collect(root) for root in collect(Iterators.product(r...))])
    return (G, roots)
  end
end
