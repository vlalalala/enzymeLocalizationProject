#%%
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
    print("y shape", dydr.shape)
    return dydr

# Boundary conditions
def bc(ya, yb):
    print("ya, yb shape", ya.shape, yb.shape)
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
    print(bc.shape)
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
# %%
membrane_radii = [1e-6, 1e-5]
R_outer = 5e-5
diffusivities = [1,1,1,1]
permeabilities = [1,1,1]
r_parts = []
r_breaks = [1e-6] + membrane_radii + [R_outer]
for i in range(len(r_breaks) - 1):
    r_segment = np.linspace(r_breaks[i], r_breaks[i+1], 100, endpoint=(i == len(r_breaks) - 2))
    if i > 0:
        r_segment = r_segment[1:]  # remove duplicate point
    r_parts.append(r_segment)

r = np.concatenate(r_parts)
r = np.unique(r)
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
# %%
import numpy as np
from scipy.integrate import solve_bvp
import matplotlib.pyplot as plt

# Physical parameters
D1 = 1.0
D2 = 0.5
P = 10.0  # Membrane permeability

r0, r1, r2 = 1.0, 2.0, 3.0
u0 = 1.0
u2 = 0.0

# ODE system: y[0] = u, y[1] = du/dr
def odes_region1(r, y):
    dydr = np.zeros_like(y)
    dydr[0] = y[1]
    dydr[1] = (-2 / r) * y[1]  # Spherical symmetry, no reaction
    return dydr

def odes_region2(r, y):
    dydr = np.zeros_like(y)
    dydr[0] = y[1]
    dydr[1] = (-2 / r) * y[1]
    return dydr

# Boundary + interface conditions
def bc(ya, yb, y_interface_left, y_interface_right):
    """
    ya: y at r = r0
    yb: y at r = r2
    y_interface_left: y at r = r1 from the left (Region 1)
    y_interface_right: y at r = r1 from the right (Region 2)
    """
    uL, duL = y_interface_left
    uR, duR = y_interface_right

    return np.array([
        ya[0] - u0,                 # u(r0) = u0
        yb[0] - u2,                 # u(r2) = u2
        D1 * duL - D2 * duR,        # Flux continuity
        -D1 * duL - P * (uR - uL)   # Jump condition
    ])

# Mesh for each region
r1_mesh = np.linspace(r0, r1, 50)
r2_mesh = np.linspace(r1, r2, 50)

# Initial guesses
y_guess_r1 = np.vstack((
    np.linspace(u0, 0.5, r1_mesh.size),
    np.zeros(r1_mesh.size)
))
y_guess_r2 = np.vstack((
    np.linspace(0.5, u2, r2_mesh.size),
    np.zeros(r2_mesh.size)
))

# Solve
sol = solve_bvp(
    fun=[odes_region1, odes_region2],
    bc=bc,
    x=[r1_mesh, r2_mesh[1:]],
    y=[y_guess_r1, y_guess_r2[1:]],
    verbose=2
)

# Plot result
if sol.success:
    r_full = np.hstack(sol.x)
    u_full = np.hstack(sol.y[0])
    plt.plot(r_full, u_full, label='u(r)')
    plt.axvline(r1, color='gray', linestyle='--', label='membrane')
    plt.xlabel('r')
    plt.ylabel('Concentration u(r)')
    plt.title('Solution with Semi-permeable Membrane')
    plt.legend()
    plt.grid(True)
    plt.show()
else:
    print("Solver failed:", sol.message)
# %%
import numpy as np
from scipy.integrate import solve_bvp
import matplotlib.pyplot as plt

# Physical parameters
D1 = 1.0
D2 = 0.5
P = 10.0  # Membrane permeability

r0, r1, r2 = 1.0, 2.0, 3.0
u0 = 1.0
u2 = 0.0

# ODE system: y[0] = u, y[1] = du/dr
def odes_region1(r, y):
    dydr = np.zeros_like(y)
    dydr[0] = y[1]
    dydr[1] = (-2 / r) * y[1]  # Spherical symmetry, no reaction
    return dydr

def odes_region2(r, y):
    dydr = np.zeros_like(y)
    dydr[0] = y[1]
    dydr[1] = (-2 / r) * y[1]
    return dydr

# Boundary + interface conditions
def bc(ya, yb, y_interface_left, y_interface_right):
    """
    ya: y at r = r0
    yb: y at r = r2
    y_interface_left: y at r = r1 from the left (Region 1)
    y_interface_right: y at r = r1 from the right (Region 2)
    """
    uL, duL = y_interface_left
    uR, duR = y_interface_right

    return np.array([
        ya[0] - u0,                 # u(r0) = u0
        yb[0] - u2,                 # u(r2) = u2
        D1 * duL - D2 * duR,        # Flux continuity
        -D1 * duL - P * (uR - uL)   # Jump condition
    ])

# Mesh for each region
r1_mesh = np.linspace(r0, r1, 50)
r2_mesh = np.linspace(r1, r2, 50)
# Initial guesses
y_guess_r1 = np.vstack((
    np.linspace(u0, 0.5, r1_mesh.size),
    np.zeros(r1_mesh.size)
))
y_guess_r2 = np.vstack((
    np.linspace(0.5, u2, r2_mesh.size),
    np.zeros(r2_mesh.size)
))
print()
#%%
r_mesh = np.append(r1_mesh, r2_mesh[1:])
r_mesh
#%%
np.ndarray.flatten(np.array([r1_mesh, r2_mesh[1:]]))
#%%
# Solve
sol = solve_bvp(
    fun=[odes_region1, odes_region2],
    bc=bc,
    x=np.append(r1_mesh, r2_mesh[1:]),
    y=np.append(y_guess_r1, y_guess_r2[1:]),
    verbose=2
)

# Plot result
if sol.success:
    r_full = np.hstack(sol.x)
    u_full = np.hstack(sol.y[0])
    plt.plot(r_full, u_full, label='u(r)')
    plt.axvline(r1, color='gray', linestyle='--', label='membrane')
    plt.xlabel('r')
    plt.ylabel('Concentration u(r)')
    plt.title('Solution with Semi-permeable Membrane')
    plt.legend()
    plt.grid(True)
    plt.show()
else:
    print("Solver failed:", sol.message)

# %%
