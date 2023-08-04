using LinearAlgebra
using TypedPolynomials
using Plots

# Define start system based on total degree
function start_system(F)
  degrees = [maxdegree(p) for p in F]
  #  @polyvar h
  #  G = [x_i^d - h^d for (d, x_i) in zip(degrees, variables(F))]
  G = [x_i^d - 1 for (d, x_i) in zip(degrees, variables(F))]
  r = [[exp(2im*pi/d)^k for k=0:d-1] for d in degrees]
  #  roots = vec([vcat(collect(root), 1) for root in collect(Iterators.product(r...))])
  roots = vec([collect(root) for root in collect(Iterators.product(r...))])
  return (G, roots)
end

function homogenize(F)
  @polyvar h
  return  [sum([h^(maxdegree(p)-maxdegree(t))*t for t in p.terms]) for p in F]
end

# Define homotopy function
function homotopy(F, G)
  γ = cis(2π * rand())
  function H(t)
    return [(1 - t) * f + γ * t * g for (f, g) in zip(F, G)]
  end
  return H
end

# Euler-Newton predictor-corrector
function en_step(H, x, t, step_size)

  # Predictor step
  vars = variables(H(t))
  # Jacobian of H evaluated at (x,t)
  JH = [jh(vars=>x) for jh in differentiate(H(t), vars)]
  # ∂H/∂t is the same as γG-F=H(1)-H(0) for our choice of homotopy
  Δx = JH \ -[gg(vars=>x) for gg in H(1)-H(0)]
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

# Adaptive step size
function adapt_step(x, x_old, step, m)
  Δ = LinearAlgebra.norm(x - x_old)
#  function adapt_step(H, x, t, step, m)
  #  Δ = LinearAlgebra.norm([h(variables(H(t))=>x) for h in H(t)])
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

# Main homotopy continuation loop
function solve(F, (G, roots) = start_system(F), maxsteps=10000)
  #  F=homogenize(F)
  H=homotopy(F,G)
  solutions = []

  for r in roots
    t = 1.0
    step_size = 0.01
    x0 = r
    m = 0

    while t > 0 && maxsteps > 0
      x = en_step(H, x0, t, step_size)
      (m, step_size) = adapt_step(x, x0, step_size, m)
      #  (m, step_size) = adapt_step(H, x, t, step_size, m)
      x0 = x
      t -= step_size
      maxsteps -= 1
    end
    push!(solutions, x0)
  end

  return solutions
end

function plot_real(solutions, F, h, v, name)
  p=plot(xlim = (-h, h), ylim = (-v, v), aspect_ratio = :equal)
  contour!(-h:0.1:h, -v:0.1:v, (x,y)->F[1](variables(F)=>[x,y]), levels=[0], cbar=false, color=:cyan)
  contour!(-h:0.1:h, -v:0.1:v, (x,y)->F[2](variables(F)=>[x,y]), levels=[0], cbar=false, color=:green)
  scatter!([real(sol[1]) for sol in solutions], [real(sol[2]) for sol in solutions], color = "red", label = "Real solutions")

  png("solutions" * name)
end

# Input polynomial system
@polyvar x y
F = [x*y - 1, x^2 + y^2 - 4]
T = [x*y - 1, x^2 + y^2 - 2]
C = [x^3 - y + 5x^2 - 10, 2x^2 - y - 10]
P = [x*y - 1, x*y]

sF = filter(u -> imag(u[1]) < 0.1 && imag(u[2]) < 0.1, solve(F))
sT = filter(u -> imag(u[1]) < 0.1 && imag(u[2]) < 0.1, solve(T))
sC = filter(u -> imag(u[1]) < 0.1 && imag(u[2]) < 0.1, solve(C))
#  sP = filter(u -> imag(u[1]) < 0.1 && imag(u[2]) < 0.1, solve(P))

# Plotting the system and the real solutions
plot_real(sF, F, 4, 4, "1")
plot_real(sT, T, 4, 4, "2")
plot_real(sC, C, 6, 12, "3")
#  plot_real(sP, P, 5, 5, "4")
