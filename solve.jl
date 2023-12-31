# External deps
using LinearAlgebra
using TypedPolynomials
using Distributed, SlurmClusterManager
slurm_manager = SlurmManager()
addprocs(slurm_manager)

# Local deps
include("random-poly.jl")
include("plot.jl")
using .RandomPoly
using .Plot
@everywhere begin
  include("start-system.jl")
  include("homotopy.jl")
  include("euler-newton.jl")
  include("adapt-step.jl")
end
# Macros defined in an @everywhere block aren't available inside it
@everywhere begin
  using .StartSystem
  using .Homotopy
  using .EulerNewton
  using .AdaptStep
end

@everywhere function compute_root(H, r, maxsteps=200)
  t = 1.0
  step_size = 0.001
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
function solve(F, G, roots)
  H = homotopy(F, G)

  result = Array{Future}(undef, length(roots))
  for i in eachindex(roots)
    result[i] = @spawnat :any compute_root(H, roots[i])
  end

  sols = Array{ComplexF64,2}(undef, length(roots), length(F))
  steps = Array{Int64}(undef, length(roots))
  for i in eachindex(roots)
    (solution, step_array) = fetch(result[i])
    sols[i, :] = solution
    steps[i] = step_array
  end

  return (sols, steps)
end

#  @polyvar x y
#  C = [x^3 - y + 5x^2 - 10, 2x^2 - y - 10]
#  Q = [x^2 + 2y, y - 3x^3]
#  F = [x*y - 1, x^2 + y^2 - 4]
#  T = [x*y - 1, x^2 + y^2 - 2]

R = random_system(3, 5)
println("System: ", R)
(G, roots)=start_system(R)
println("Number of roots: ", length(roots))

# Parallel execution
println("PARALLEL")
@time begin
  (sol, steps) = solve(R, G, roots)
end
println("Number of steps: ", steps)
# converting sR to array of arrays instead of a matrix
sol = [sol[i, :] for i in 1:length(sol[:, 1])]
sol = filter(u -> imag(u[1]) < 0.1 && imag(u[2]) < 0.1, sol)
sol = map(u -> real.(u), sol)
vars = variables(R)
println("Solutions: ", sol)
println("Norms (lower = better): ", [norm([f(vars => s) for f in R]) for s in sol])

# Single execution
println("SINGLE")
wait(rmprocs(workers()))
@time begin
  (sol, steps) = solve(R, G, roots)
end
println("Number of steps: ", steps)
# converting sR to array of arrays instead of a matrix
sol = [sol[i, :] for i in 1:length(sol[:, 1])]
sol = filter(u -> imag(u[1]) < 0.1 && imag(u[2]) < 0.1, sol)
sol = map(u -> real.(u), sol)
vars = variables(R)
println("Solutions: ", sol)
println("Norms (lower = better): ", [norm([f(vars => s) for f in R]) for s in sol])

# See https://github.com/kleinhenz/SlurmClusterManager.jl/issues/11
finalize(slurm_manager)

# Plotting the system and the real solutions
# ENV["GKSwstype"] = "nul"
# plot_real(sC, C, 6, 12, "1")
# plot_real(sQ, Q, 2, 2, "2")
# plot_real(sF, F, 4, 4, "3")
# plot_real(sT, T, 4, 4, "4")
# plot_real(sol, R, 5, 5, "random")
