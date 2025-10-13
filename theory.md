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

and at $r=R$ we have

$$ \frac{du}{dr}(R) = [u_\mathrm{ext} − u(R)] \cdot \frac{p_\mathrm{u}}{D_\mathrm{u}}$$

with $p_\mathrm{u}$ the permeability to the outer membrane.

Flux through the membrane using the permeability law

$$J_\mathrm{membrane} = p \cdot (u_R - u_L)

## With multiple regions

In each region $k$ we have 

$$ \frac{1}{r^2} \frac{d}{dr} \left( r^2 \frac{du_k}{dr} \right) + \frac{1}{D_k} R(u) = 0 ,  $$

$$ D_k \, \frac{\mathrm{d}u_k}{\mathrm{d}r} \big|_{r = R_k^-}
= D_{k+1} \, \frac{\mathrm{d}u_{k+1}}{\mathrm{d}r} \big|_{r = R_k^+}
= p_k \left( u_{k+1}(R_k^+) - u_k(R_k^-) \right)$$

```python
import numpy as np
from scipy.integrate import solve_bvp
import matplotlib.pyplot as plt

# --- Parameters ---
R1 = 1.0     # membrane radius
R2 = 2.0     # outer radius
D1 = 1.0     # diffusion in region 1
D2 = 0.5     # diffusion in region 2
p = 0.3      # permeability at membrane
u_ext = 1.0  # concentration outside

# --- Define mesh ---
r1 = np.linspace(1e-6, R1, 100)        # avoid r = 0
r2 = np.linspace(R1, R2, 100)
r = np.concatenate((r1, r2))
regions = np.concatenate((np.zeros_like(r1), np.ones_like(r2)))  # 0 for inner, 1 for outer

# --- Define ODE system ---
def fun(r, y):
    dydr = np.zeros_like(y)
    dydr[0] = y[1]                            # y = u, z = du/dr
    dydr[1] = -2 * y[1] / r                   # from spherical coordinates
    return dydr

# --- Apply different D in regions ---
def wrapped_fun(r, y):
    dydr = fun(r, y)
    D = np.where(r <= R1, D1, D2)
    dydr[0] /= D
    dydr[1] /= D
    return dydr

# --- Boundary conditions ---
def bc(ya, yb):
    bc_vals = np.zeros(4)
    # ya: values at r = r[0] ~ 0 (inner)
    # yb: values at r = r[-1] = R2 (outer)
    # inner symmetry, outer boundary, flux continuity at the membrane from the left, flux continuity at the membrane from the right
    
    # At r=0: symmetry → du/dr = 0
    bc_vals[0] = ya[1] # du/dr (at r=0) = 0
    
    # At r=R2: Dirichlet condition
    bc_vals[1] = yb[0] - u_ext

    # Find index of the mesh point closest to membrane R1
    idx = np.where(np.isclose(r, R1))[0][0]
    
    # Interface condition: flux continuity
    u_L, du_L = sol.sol(R1 - 1e-6)  # left limit
    u_R, du_R = sol.sol(R1 + 1e-6)  # right limit
    # Gets the values of the solution just to the left and right of the membrane
    # Flux through the membrane using the permeability law
    J_membrane = p * (u_R - u_L)
    
    # Flux balance at membrane
    bc_vals[2] = D1 * du_L - J_membrane
    bc_vals[3] = D2 * du_R - J_membrane
    # Enforce that fluxes on both sides of the membrane equal the membrane flux

    return bc_vals

# --- Initial guess: linear profile ---
y_init = np.zeros((2, r.size))
y_init[0] = np.linspace(0, u_ext, r.size)

# --- Solve ---
sol = solve_bvp(wrapped_fun, bc, r, y_init)

# --- Plot solution ---
if sol.success:
    plt.plot(sol.x, sol.y[0], label="u(r)")
    plt.xlabel("r")
    plt.ylabel("Concentration u")
    plt.title("Spherical diffusion with a membrane at r = R1")
    plt.grid()
    plt.legend()
    plt.show()
else:
    print("Solution failed:", sol.message)
```

Correct:?:
```python
import numpy as np
from scipy.integrate import solve_bvp
import matplotlib.pyplot as plt

# Parameters
R1 = 1.0    # membrane position
R2 = 2.0    # outer boundary
D1 = 1.0    # diffusion coefficient inner region
D2 = 0.5    # diffusion coefficient outer region
p = 0.3     # membrane permeability
u_ext = 1.0 # external concentration

# Mesh: ensure R1 is in the mesh
r1 = np.linspace(1e-6, R1, 100)
r2 = np.linspace(R1, R2, 100)[1:]  # avoid duplicate R1
r = np.concatenate((r1, r2))

# Define full system: y = [u1, u1', u2, u2']
def fun(r, y):
    dydr = np.zeros_like(y)
    dydr[0] = y[1]
    dydr[1] = -2 * y[1] / r
    dydr[2] = y[3]
    dydr[3] = -2 * y[3] / r
    return dydr

# Boundary conditions
def bc(ya, yb):
    bc = np.zeros(4)
    
    # At r=0: symmetry
    bc[0] = ya[1]  # u1'(0) = 0
    
    # At r=R2: Dirichlet
    bc[1] = yb[2] - u_ext  # u2(R2) = u_ext

    # At r=R1: interface conditions
    u1_R1 = ya[0]   # u1 at R1
    u1p_R1 = ya[1]  # du1/dr at R1
    u2_R1 = yb[2]   # u2 at R1
    u2p_R1 = yb[3]  # du2/dr at R1

    J_membrane = p * (u2_R1 - u1_R1)
    bc[2] = D1 * u1p_R1 - J_membrane
    bc[3] = D2 * u2p_R1 - J_membrane

    return bc

# Initial guess: linear from 0 to u_ext
y_guess = np.zeros((4, r.size))
y_guess[0, :] = np.linspace(0, u_ext, r.size)  # u1 guess
y_guess[2, :] = np.linspace(0, u_ext, r.size)  # u2 guess

# Solve
sol = solve_bvp(fun, bc, r, y_guess)

# Plot result
if sol.success:
    plt.plot(sol.x, sol.y[0], label='u1 (inner region)')
    plt.plot(sol.x, sol.y[2], label='u2 (outer region)')
    plt.axvline(R1, color='gray', linestyle='--', label='Membrane at r=R1')
    plt.xlabel('r')
    plt.ylabel('u(r)')
    plt.legend()
    plt.grid()
    plt.title('Steady-state diffusion with membrane at r=R1')
    plt.show()
else:
    print("Solver failed:", sol.message)

```
So for 2 membranes (3 regions), your state is:
```python

y = [u0, u0', u1, u1', u2, u2']
```
For each region 
k∈[0,N]
k∈[0,N], track:

concentration u_k and radial derivative u_k'.

```
def fun(r, y):
    dydr = np.zeros_like(y)
    for k in range(num_regions):
        u = y[2 * k]
        up = y[2 * k + 1]
        dydr[2 * k] = up
        dydr[2 * k + 1] = -2 * up / r  # spherical term
    return dydr
```

bc(ya, yb)
first condition u'_0(0) = 0, last condition u_N(R_N+1) = u_ext

For each membrane  add:

Flux jump condition:

```
Dk * u_k'(R_k^-) - J_membrane = 0
Dk+1 * u_{k+1}'(R_k^+) - J_membrane = 0
```

FLUX:
$$J = p_k​(u_{k+1}​(R_{k+​})−u_k​(R_{k−​}))$$


Complete
```python
import numpy as np
from scipy.integrate import solve_bvp
import matplotlib.pyplot as plt

# --- Parameters ---
membrane_radii = [1.0, 2.0]          # Membrane positions: R1, R2
R_outer = 3.0                        # Outer boundary
diffusivities = [1.0, 0.5, 0.2]      # D0, D1, D2 for regions 0, 1, 2
permeabilities = [0.3, 0.1]          # p0 (R1), p1 (R2)
u_ext = 1.0                          # External concentration at r=R_outer

# --- Mesh setup ---
r_parts = []
r_breaks = [1e-6] + membrane_radii + [R_outer]
for i in range(len(r_breaks) - 1):
    r_segment = np.linspace(r_breaks[i], r_breaks[i+1], 100, endpoint=(i == len(r_breaks) - 2))
    if i > 0:
        r_segment = r_segment[1:]  # remove duplicate point
    r_parts.append(r_segment)

r = np.concatenate(r_parts)
num_regions = len(diffusivities)

# --- ODE system ---
def fun(r, y):
    dydr = np.zeros_like(y)
    for k in range(num_regions):
        u = y[2*k]
        up = y[2*k + 1]
        dydr[2*k] = up
        dydr[2*k + 1] = -2 * up / r
    return dydr

# --- Boundary and Interface Conditions ---
def bc(ya, yb):
    bc_vals = []

    # 1. Symmetry at r = 0
    bc_vals.append(ya[1])  # u0'(0) = 0

    # 2. Dirichlet at outer boundary
    bc_vals.append(yb[-2] - u_ext)  # uN(R_outer) = u_ext

    # 3. Interface conditions at each membrane
    for i, Rm in enumerate(membrane_radii):
        # Left side values (region i)
        u_L, du_L = sol.sol(Rm - 1e-6)[2*i:2*i+2]
        # Right side values (region i+1)
        u_R, du_R = sol.sol(Rm + 1e-6)[2*(i+1):2*(i+1)+2]

        J = permeabilities[i] * (u_R - u_L)

        bc_vals.append(diffusivities[i] * du_L - J)
        bc_vals.append(diffusivities[i+1] * du_R - J)

    return np.array(bc_vals)

# --- Initial guess ---
y_guess = np.zeros((2 * num_regions, r.size))
for k in range(num_regions):
    y_guess[2*k] = np.linspace(0, u_ext, r.size)

# --- Solve ---
sol = solve_bvp(fun, bc, r, y_guess)

# --- Plot results ---
if sol.success:
    plt.figure(figsize=(8, 5))
    for k in range(num_regions):
        plt.plot(sol.x, sol.y[2*k], label=f'u_{k}(r)')
    for R in membrane_radii:
        plt.axvline(R, color='gray', linestyle='--', label=f'Membrane at r={R}')
    plt.xlabel('r')
    plt.ylabel('u(r)')
    plt.title('Steady-State Diffusion with Multiple Membranes')
    plt.grid()
    plt.legend()
    plt.tight_layout()
    plt.show()
else:
    print("Solver failed:", sol.message)
```

# Typical values

Diffusion coefficient of l-Trp taken as $6.6\cdot 10^{-6} cm^{2} s^{-1} = 6.6\cdot 10^{-10} m^{2} s^{-1}$ (estimated from https://pmc.ncbi.nlm.nih.gov/articles/PMC16526/).
Assuming same diffusion coefficient for the other values.

Following https://pubs.acs.org/doi/pdf/10.1021/acsomega.3c08233?ref=article_openPDF (quantification of kinetics of VioA) we assume $k_\mathrm{M} = 125 \mu M = 1.25\cdot 10^{-4}M$ and $k_\mathrm{cat} = 0.75 s^{-1}$ for all enzymes.

For the external concentration of l-Trp we're assuming $s_0 = 25 \mu M = 25\cdot 10^{-6}M$ and the permeability $p = 90 \mu M s^{-1} = 90\cdot 10^{-6}M s^{-1}$.
(as in https://www.cell.com/biophysj/fulltext/S0006-3495(16)34263-1?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0006349516342631%3Fshowall%3Dtrue)

We use a radius of $1\mu m$.

The concentration of enzymes is chosen as $25 mM = 25 \cdot 10^{-3}M$

To find the total number of moles $n$ within the volume with radius R, calculate
$$
n = \int_0^R C(r) \cdot 4\pi r^2 \mathrm{d}r
$$
Important: convert $C(r)$ from $M$ (moles per litre) to moles per $m^{3}$: $M = \frac{moles}{m^3}\cdot 10^{-3}$.

So 
$$
n = \int_0^R C(r) \cdot 10^{3}\cdot 4\pi r^2 \mathrm{d}r
$$
with n in moles and C(r) in M.
