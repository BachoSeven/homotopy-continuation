# External dependencies
using TypedPolynomials

# Local dependencies
include("homotopy.jl")
include("plot.jl")
include("euler-newton.jl")
include("adapt-step.jl")
include("start-system.jl")
include("homogenize.jl")
using .Homotopy
using .Plot
using .EulerNewton
using .AdaptStep
using .StartSystem
using .Homogenize

# Main homotopy continuation loop
function solve(F, (G, roots) = start_system(F), maxsteps=10000)
  #  F=homogenize(F)
  H=homotopy(F,G)
  solutions = []
  steps = 0

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

  return (solutions, steps)
end

# Input polynomial system
@polyvar x y
F = [x*y - 1, x^2 + y^2 - 4]
T = [x*y - 1, x^2 + y^2 - 2]
C = [x^3 - y + 5x^2 - 10, 2x^2 - y - 10]

(sF, sf) = solve(F)
(sT, st) = solve(T)
(sC, sc) = solve(C)

println(sf)
println(st)
println(sc)

sF = filter(u -> imag(u[1]) < 0.1 && imag(u[2]) < 0.1, sF)
sT = filter(u -> imag(u[1]) < 0.1 && imag(u[2]) < 0.1, sT)
sC = filter(u -> imag(u[1]) < 0.1 && imag(u[2]) < 0.1, sC)

# Plotting the system and the real solutions
ENV["GKSwstype"]="nul"
plot_real(sF, F, 4, 4, "1")
plot_real(sT, T, 4, 4, "2")
plot_real(sC, C, 6, 12, "3")
