# External dependencies
using TypedPolynomials
using LinearAlgebra
using Distributed
using ClusterManagers

# Local dependencies
include("random_poly.jl")
include("start-system.jl")
include("homotopy.jl")
include("euler-newton.jl")
include("adapt-step.jl")
include("plot.jl")
using .RandomPoly
using .StartSystem
using .Homotopy
using .EulerNewton
using .AdaptStep
using .Plot

# Launch worker processes
addprocs(SlurmManager(40), N=20, t="01:00:00"))

function compute_root(H, r, maxsteps)
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
  return (x0, steps)
end

# Main homotopy continuation loop
function solve(F, (G, roots) = start_system(F))
  H=homotopy(F,G)

  @distributed for r in roots
    (solutions, step_array) = compute_root(H, r, maxsteps = 1000)
  end

  # Gather results from worker processes
  sols = fetch(solutions)
  steps = fetch(step_array)
  return (sols, steps)
end

# Input polynomial systems
#  @polyvar x y
#  C = [x^3 - y + 5x^2 - 10, 2x^2 - y - 10]
#  Q = [x^2 + 2y, y - 3x^3]
#  F = [x*y - 1, x^2 + y^2 - 4]
#  T = [x*y - 1, x^2 + y^2 - 2]
dimension = 2
R = random_system(2, 2)
println(R)

#  (sC, stepsC) = solve(C)
#  (sQ, stepsQ) = solve(Q)
#  (sF, stepsF) = solve(F)
#  (sT, stepsT) = solve(T)
(sR, stepsR) = solve(R)

#  println("C: ", stepsC)
#  println("Q: ", stepsQ)
#  println("F: ", stepsF)
#  println("T: ", stepsT)
println("R: ", stepsR)

#  sC = filter(u -> imag(u[1]) < 0.1 && imag(u[2]) < 0.1, sC)
#  sQ = filter(u -> imag(u[1]) < 0.1 && imag(u[2]) < 0.1, sQ)
#  sF = filter(u -> imag(u[1]) < 0.1 && imag(u[2]) < 0.1, sF)
#  sT = filter(u -> imag(u[1]) < 0.1 && imag(u[2]) < 0.1, sT)
sR = filter(u -> imag(u[1]) < 0.1 && imag(u[2]) < 0.1, sR)

vars = variables(R)
println("solutions: ", sR)
println([LinearAlgebra.norm([f(vars=>s) for f in R]) for s in sR])

# Plotting the system and the real solutions
ENV["GKSwstype"]="nul"
#  plot_real(sC, C, 6, 12, "1")
#  plot_real(sQ, Q, 2, 2, "2")
#  plot_real(sF, F, 4, 4, "3")
#  plot_real(sT, T, 4, 4, "4")
plot_real(sR, R, 5, 5, "random")
