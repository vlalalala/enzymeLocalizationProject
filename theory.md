# Reaction-diffusion system in spherical coordinates

In spherical symmetry (assuming no angular dependence), the steady-state reaction-diffusion equation becomes

$$ \frac{1}{r^2} \frac{d}{dr} \left( r^2 \frac{du}{dr} \right) + \frac{1}{D} R(u) = 0 ,  $$

which simplifies to

$$\frac{d^2u}{dr^2} + \frac{2}{r} \frac{du}{dr} + \frac{1}{D} R(u) = 0$$

(Singularity at r=0!)

Using
$$\begin{cases}
y_0 = u \\
y_1 = \frac{du}{dr} \\\end{cases}
$$

this leads to 
$$
\begin{cases}
\frac{dy_0}{dr} = y_1 \\
\frac{dy_1}{dr} = -\frac{2}{r} y_1 - \frac{1}{D} R(y_0)
\end{cases}$$

As for the boundary conditions:

At $r=0$ we use reflection

$$\frac{du}{dr}(0) = 0$$

and at $r=R$ we set a Dirichlet condition

$$u(R) = u_R $$

```python
import numpy as np
from scipy.integrate import solve_bvp
import matplotlib.pyplot as plt

# Parameters
D = 1.0  # diffusion coefficient

# Reaction term R(u), define as needed
def reaction(u):
    return u * (1 - u)  # example: logistic growth

# First-order system of ODEs
def odes(r, y):
    y0, y1 = y
    # Avoid division by zero at r = 0
    with np.errstate(divide='ignore', invalid='ignore'):
        dy1 = -2 / r * y1 - reaction(y0) / D
        dy1 = np.where(r == 0, -reaction(y0) / D, dy1)
    return np.vstack((y1, dy1))

# Boundary conditions: u(0) finite (du/dr = 0), u(R) = u_R
def bc(ya, yb):
    return np.array([ya[1], yb[0] - u_R])  # ya[1]=du/dr at r=0, yb[0]=u(R)

# Domain
R = 1.0  # outer radius
r = np.linspace(0, R, 100)

# Boundary value at outer radius
u_R = 1.0

# Initial guess
y_guess = np.zeros((2, r.size))
y_guess[0] = u_R * (1 - (r / R)**2)  # trial function satisfying u(R) = u_R

# Solve BVP
sol = solve_bvp(odes, bc, r, y_guess)

# Plot solution
if sol.success:
    r_plot = np.linspace(0, R, 200)
    u_plot = sol.sol(r_plot)[0]
    plt.plot(r_plot, u_plot, label='u(r)')
    plt.xlabel('r')
    plt.ylabel('u')
    plt.title('Steady-State Reaction-Diffusion in Sphere')
    plt.grid(True)
    plt.legend()
    plt.show()
else:
    print("BVP solution failed.")

```