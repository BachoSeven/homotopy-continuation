module RandomPoly
  export random_system

  using TypedPolynomials
  using Random
  using Distributions

  # Random polynomial of degree n in m variables
  function random_poly(n, m)
    x = [TypedPolynomials.Variable{Symbol("x[$i]")}() for i in 1:m]

    monomial_powers=collect(Iterators.product([0:n for _ in 1:m]...))
    monomials = [prod(x.^i) for i in monomial_powers if sum(i) == n]

    return sum(map(m -> rand(Normal()) * m, monomials))
  end

  # Generate a system of m random polynomials in m variables
  # of degree d_i randomly chosen between 1 and max_degree
  function random_system(m, max_degree)
    d = rand(1:max_degree, m)
    println("Expected number of real zeros: ", sqrt(prod(d)))
    random_polys = [random_poly(d[i], m) for i in 1:m]

    return random_polys
  end
end
