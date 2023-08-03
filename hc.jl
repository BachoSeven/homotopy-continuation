using LinearAlgebra
using TypedPolynomials
using Plots

# Define start system based on total degree
function start_system(F)
  degrees = [maxdegree(p) for p in F]
  G = [x_i^d - 1 for (d, x_i) in zip(degrees, variables(F))]
  r = [[exp(2im*pi/d)^k for k=0:d-1] for d in degrees]
  roots = collect(Iterators.product(r...))
  return (G, roots)
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
  for _ in 1:5
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
function adapt_step(H, x, t, step, m)
  Δ = LinearAlgebra.norm([h(variables(H(t))=>x) for h in H(t)])
  if Δ > 0.1
    step = 0.5 * step
  elseif Δ < 0.001
    m+=1
    if (m == 5)
      step = 2 * step
      m = 0
    end
  end

  return (m, step)
end

# Main homotopy continuation loop
function solve(F, maxsteps=10000)
  (G, roots) = start_system(F)
  H=homotopy(F,G)
  solutions = []

  for r in roots
    t = 1.0
    step_size = 0.1
    x0 = r
    m = 0

    while t > 0 && maxsteps > 0
      x0 = en_step(H, x0, t, step_size)
      (m, step_size) = adapt_step(H, x0, t, step_size, m)
      t -= step_size
      maxsteps -= 1
    end
    push!(solutions, x0)
  end

  return solutions
end

function plot_real(solutions, F)
  p=plot(xlim = (-3, 3), ylim = (-3, 3), aspect_ratio = :equal)
  contour!(-3:0.1:3, -3:0.1:3, (x,y)->F[1](variables(F)=>[x,y]), levels=[0], cbar=false, color=:cyan)
  contour!(-3:0.1:3, -3:0.1:3, (x,y)->F[2](variables(F)=>[x,y]), levels=[0], cbar=false, color=:green)
  scatter!([real(sol[1]) for sol in solutions], [real(sol[2]) for sol in solutions], color = "red", label = "Solutions")

  png("solutions")
end

# Input polynomial system
@polyvar x y
F = [x*y - 1, x^2 + y^2 - 4]

sF = solve(F)

# Plotting the system and the real solutions
plot_real(sF, F)
