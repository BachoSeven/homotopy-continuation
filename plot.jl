module Plot
  using Plots, TypedPolynomials

  export plot_real

  function plot_real(solutions, F, h, v, name)
    plot(xlim = (-h, h), ylim = (-v, v), aspect_ratio = :equal)
    contour!(-h:0.1:h, -v:0.1:v, (x,y)->F[1](variables(F)=>[x,y]), levels=[0], cbar=false, color=:cyan)
    contour!(-h:0.1:h, -v:0.1:v, (x,y)->F[2](variables(F)=>[x,y]), levels=[0], cbar=false, color=:green)
    scatter!([real(sol[1]) for sol in solutions], [real(sol[2]) for sol in solutions], color = "red", label = "Real solutions")

    png(joinpath("./plots", "solutions" * name))
  end
end
