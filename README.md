# Homotopy Continuation in Julia

This is a project for the "Laboratorio Computazionale" exam at the University of Pisa

## Implemented

- Total-degree Homotopy with "Roots of unity" start system
- Euler-Newton predictor-corrector method with adaptive step size
- Homotopy Continuation for all roots of the target system

## TODO

- Projective coordinates
- Parallelization
- Extract functions in separate modules(?)

## Example systems

Here's some tests on 2x2 systems, with the plotted real approximate solutions

$$
\begin{align*}
x^2 + y^2 - 4 &= 0 \\
xy - 1 &= 0 \\
\end{align*}
$$

![](solutions1.png)

---

$$
\begin{align*}
x^2 + y^2 - 2 &= 0 \\
xy - 1 &= 0 \\
\end{align*}
$$

![](solutions2.png)

---

$$
\begin{align*}
x^3 + 5x^2 - y - 10 &= 0 \\
2x^2 - y - 10 &= 0 \\
\end{align*}
$$

![](solutions3.png)
