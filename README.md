# Solving the Stokes equation using the MAC scheme and multigrid methods
## Final project @ Numerical Linear Algebra (数值代数) at Peking University. 
Consider the Stokes equations:
$$-\Delta \vec{u} + \nabla p  = \vec{F},  (x, y) \in (0,1) \times (0,1),$$
$$\text{div} \vec{u}  = 0,  (x, y) \in (0,1) \times (0,1).$$

The boundary conditions are:
$$\frac{\partial u}{\partial \vec{n}} = b, \quad y = 0, \quad \frac{\partial u}{\partial \vec{n}} = t, \quad y = 1,$$
$$\frac{\partial v}{\partial \vec{n}} = l, \quad x = 0, \quad \frac{\partial v}{\partial \vec{n}} = r, \quad x = 1,$$
$$u = 0, \quad x = 0, 1, \quad v = 0, \quad y = 0, 1.$$

Here, $\vec{u} = (u, v)$ represents the velocity, $p$ the pressure, $\vec{F} = (f, g)$ the external force, and $\vec{n}$ the outward normal direction.
In the domain $\Omega = (0,1) \times (0,1)$, the external force is:
$$f(x, y) = -4 \pi^2 (2 \cos(2 \pi x) - 1) \sin(2 \pi y) + x^2, \\
g(x, y) = 4 \pi^2 (2 \cos(2 \pi y) - 1) \sin(2 \pi x).$$

The true solution of the Stokes equations is:
$$u(x, y) = (1 - \cos(2 \pi x)) \sin(2 \pi y), \\
v(x, y) = -(1 - \cos(2 \pi y)) \sin(2 \pi x), \\
p(x, y) = \frac{x^3}{3} - \frac{1}{12}.$$

We will design a MAC scheme for $u$, $v$, and $p$ on a staggered grid and discretize the original problem using the finite difference method. Approximate the operators directly at internal grid points using the difference formula, discretize at Neumann boundaries using ghost cells, and no equations need to be formulated at enforced boundary conditions. The linear systems are solved by multigrid methods. See the report for more details (the report is written in Chinese).


