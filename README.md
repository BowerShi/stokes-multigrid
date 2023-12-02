# Solving the Stokes equation using the MAC scheme and multigrid methods
## Final project @ Numerical Linear Algebra (数值代数) at Peking University. 
Consider the Stokes equations:
$$
\begin{cases}
-\Delta \vec{u} + \nabla p & = \vec{F}, & (x, y) \in (0,1) \times (0,1), \\
\operatorname{div} \vec{u} & = 0, & (x, y) \in (0,1) \times (0,1).
\end{cases}
$$

The boundary conditions are:
$$
\begin{aligned}
&\frac{\partial u}{\partial \vec{n}} = b, \quad y = 0, \quad \frac{\partial u}{\partial \vec{n}} = t, \quad y = 1, \\
&\frac{\partial v}{\partial \vec{n}} = l, \quad x = 0, \quad \frac{\partial v}{\partial \vec{n}} = r, \quad x = 1, \\
&u = 0, \quad x = 0, 1, \quad v = 0, \quad y = 0, 1.
\end{aligned}
$$

Here, \(\vec{u} = (u, v)\) represents the velocity, \(p\) the pressure, \(\vec{F} = (f, g)\) the external force, and \(\vec{n}\) the outward normal direction.
In the domain \(\Omega = (0,1) \times (0,1)\), the external force is:
$$
\begin{aligned}
& f(x, y) = -4 \pi^2 (2 \cos(2 \pi x) - 1) \sin(2 \pi y) + x^2, \\
& g(x, y) = 4 \pi^2 (2 \cos(2 \pi y) - 1) \sin(2 \pi x).
\end{aligned}
$$

The true solution of the Stokes equations is:
$$
\begin{aligned}
& u(x, y) = (1 - \cos(2 \pi x)) \sin(2 \pi y), \\
& v(x, y) = -(1 - \cos(2 \pi y)) \sin(2 \pi x), \\
& p(x, y) = \frac{x^3}{3} - \frac{1}{12}.
\end{aligned}
$$

#### Staggered Grid MAC Scheme
Design a MAC scheme for \(u\), \(v\), and \(p\) on a staggered grid and discretize the original problem using the finite difference method. Approximate the operators directly at internal grid points using the difference formula, discretize at Neumann boundaries using ghost cells, and no equations need to be formulated at enforced boundary conditions.





Acknoledgements: Some pictures in the report are created by Prof. Jun Hu (instructor) and Qingyu Wu (TA).
