# External dependencies
using TypedPolynomials

# Local dependencies
include("start-system.jl")
include("homotopy.jl")
#  include("homogenize.jl")
include("euler-newton.jl")
include("adapt-step.jl")
include("plot.jl")
using .StartSystem
using .Homotopy
#  using .Homogenize
using .EulerNewton
using .AdaptStep
using .Plot

# Main homotopy continuation loop
function solve(F, (G, roots) = start_system(F), maxsteps = 1000)
  #  F=homogenize(F)
  H=homotopy(F,G)
  solutions = []
  step_array = []

  Threads.@threads for r in roots
    t = 1.0
    step_size = 0.01
    x0 = r
    m = 0
    steps = 0

    while t > 0 && steps < maxsteps
      x0 = en_step(H, x0, t, step_size)
      (m, step_size) = adapt_step(H, x0, t, step_size, m)
      t -= step_size
      steps += 1
    end
    push!(solutions, x0)
    push!(step_array, steps)
  end

  return (solutions, step_array)
end

# Input polynomial system
@polyvar x y
C = [x^3 - y + 5x^2 - 10, 2x^2 - y - 10]
Q = [x^2 + 2y, y - 3x^3]
F = [x*y - 1, x^2 + y^2 - 4]
T = [x*y - 1, x^2 + y^2 - 2]

(sC, stepsC) = solve(C)
(sQ, stepsQ) = solve(Q)
(sF, stepsF) = solve(F)
(sT, stepsT) = solve(T)

println("C: ", stepsC)
println("Q: ", stepsQ)
println("F: ", stepsF)
println("T: ", stepsT)

sC = filter(u -> imag(u[1]) < 0.1 && imag(u[2]) < 0.1, sC)
sQ = filter(u -> imag(u[1]) < 0.1 && imag(u[2]) < 0.1, sQ)
sF = filter(u -> imag(u[1]) < 0.1 && imag(u[2]) < 0.1, sF)
sT = filter(u -> imag(u[1]) < 0.1 && imag(u[2]) < 0.1, sT)

# Plotting the system and the real solutions
ENV["GKSwstype"]="nul"
plot_real(sC, C, 6, 12, "1")
plot_real(sQ, Q, 2, 2, "2")
plot_real(sF, F, 4, 4, "3")
plot_real(sT, T, 4, 4, "4")
