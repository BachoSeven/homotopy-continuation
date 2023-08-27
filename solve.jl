# External dependencies
using TypedPolynomials

# Local dependencies
include("start-system.jl")
include("homotopy.jl")
include("homogenize.jl")
include("euler-newton.jl")
include("adapt-step.jl")
include("plot.jl")
using .StartSystem
using .Homotopy
using .Homogenize
using .EulerNewton
using .AdaptStep
using .Plot

# Main homotopy continuation loop
function solve(F, (G, roots) = start_system(F), maxsteps=10000)
  #  F=homogenize(F)
  H=homotopy(F,G)
  solutions = []

  @time Threads.@threads for r in roots
    t = 1.0
    step_size = 0.01
    x0 = r
    m = 0
    steps = 0

    while t > 0 && steps < maxsteps
      x = en_step(H, x0, t, step_size)
      (m, step_size) = adapt_step(x, x0, step_size, m)
      x0 = x
      t -= step_size
      steps += 1
    end
    push!(solutions, x0)
  end

  return solutions
end

# Input polynomial system
@polyvar x y
F = [x*y - 1, x^2 + y^2 - 4]
T = [x*y - 1, x^2 + y^2 - 2]
C = [x^3 - y + 5x^2 - 10, 2x^2 - y - 10]

sF = filter(u -> imag(u[1]) < 0.1 && imag(u[2]) < 0.1, solve(F))
sT = filter(u -> imag(u[1]) < 0.1 && imag(u[2]) < 0.1, solve(T))
sC = filter(u -> imag(u[1]) < 0.1 && imag(u[2]) < 0.1, solve(C))

# Plotting the system and the real solutions
ENV["GKSwstype"]="nul"
plot_real(sF, F, 4, 4, "1")
plot_real(sT, T, 4, 4, "2")
plot_real(sC, C, 6, 12, "3")
