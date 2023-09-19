module RandomPoly
  export random_system

  using TypedPolynomials
  using Random
  using Distributions

  # Random polynomial of degree n in m variables
  function random_poly(n, m)
    x = [TypedPolynomials.Variable{Symbol("x[$i]")}() for i in 1:m]

    monomial_powers=vcat(collect(Iterators.product([0:n for _ in 1:m]...))...)

    return sum(map(i -> rand(Uniform(-10,10)) * prod(x.^i), monomial_powers))
  end

  # Generate a system of m random polynomials in m variables of degree d_i
  function random_system(m, max_degree)
    d = rand(1:max_degree, m)
    println("generating system")
    random_polys = [random_poly(d[i], m) for i in 1:m]
    println("done generating system")

    return random_polys
  end
end
