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
import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve

# Domain and discretization
a, b = 0, 1
N = 101
x = np.linspace(a, b, N)
h = x[1] - x[0]

# Internal condition location
c = 0.5
ic_idx = np.argmin(np.abs(x - c))  # nearest grid point

# Finite-difference matrix for y'' = -y
A = diags([1, -2, 1], [-1, 0, 1], shape=(N, N)).tocsc() / h**2
A -= diags([1], [0])  # subtract y term (i.e. y'' + y = 0)

# Right-hand side
bvec = np.zeros(N)

# Boundary conditions
A[0, :] = 0
A[0, 0] = 1
bvec[0] = 0

A[-1, :] = 0
A[-1, -1] = 1
bvec[-1] = 1

# Internal derivative jump at x = c
# Approx y'(c) ≈ (y_{i+1} - y_{i-1}) / (2h)
A[ic_idx, :] = 0
A[ic_idx, ic_idx-1] = -1/(2*h)
A[ic_idx, ic_idx+1] =  1/(2*h)
bvec[ic_idx] = 2  # jump = 2

# Solve linear system
y = spsolve(A, bvec)

# (Optional) visualize
import matplotlib.pyplot as plt
plt.plot(x, y, label="solution")
plt.axvline(c, color='r', linestyle='--', label='internal point')
plt.legend()
plt.show()
# %%
# skeleton_nlbvp_membrane.py
import numpy as np
import dedalus.public as d3

# --- parameters (example) ---
Nr = 64               # resolution per patch
r0, r1, r2, r3 = 0.0, 0.5, 1.0, 2.0   # radii: membranes at r1 and r2 etc.
D1, D2, D3 = 1.0, 0.5, 1.2             # diffusion in each region
k12, k23 = 1.0, 0.2                    # membrane conductances at r=r1 and r=r2
dealias = 2
dtype = np.float64
tolerance = 1e-10

# --- Coordinates / domains (1D radial patches with Chebyshev) ---
r1coord = d3.Coordinate('r1')
dist1   = d3.Distributor(r1coord, dtype=dtype)
basis1  = d3.Chebyshev(r1coord, Nr, bounds=(r0, r1), dealias=dealias)

r2coord = d3.Coordinate('r2')
dist2   = d3.Distributor(r2coord, dtype=dtype)
basis2  = d3.Chebyshev(r2coord, Nr, bounds=(r1, r2), dealias=dealias)

r3coord = d3.Coordinate('r3')
dist3   = d3.Distributor(r3coord, dtype=dtype)
basis3  = d3.Chebyshev(r3coord, Nr, bounds=(r2, r3), dealias=dealias)

# --- Fields (unknown concentrations on each patch) ---
c1 = dist1.Field(name='c1', bases=basis1)
c2 = dist2.Field(name='c2', bases=basis2)
c3 = dist3.Field(name='c3', bases=basis3)

# tau fields for BC enforcement (tau basis degree depends on highest derivative)
tau1 = dist1.Field(name='tau1', bases=basis1.derivative_basis(2))
tau2 = dist2.Field(name='tau2', bases=basis2.derivative_basis(2))
tau3 = dist3.Field(name='tau3', bases=basis3.derivative_basis(2))

# --- helpful substitutions/operators ---
dr1 = lambda A: d3.Differentiate(A, r1coord)
dr2 = lambda A: d3.Differentiate(A, r2coord)
dr3 = lambda A: d3.Differentiate(A, r3coord)

# radial laplacian in spherical symmetry: d2/dr2 + (2/r) d/dr
lap1 = lambda A: d3.Differentiate(A, r1coord, 2) + (2/d3.Coordinate('r1')) * dr1(A)
lap2 = lambda A: d3.Differentiate(A, r2coord, 2) + (2/d3.Coordinate('r2')) * dr2(A)
lap3 = lambda A: d3.Differentiate(A, r3coord, 2) + (2/d3.Coordinate('r3')) * dr3(A)

# Reaction terms (example nonlinear functions)
R1 = lambda A: 0.0*A         # replace with actual reaction, e.g. k*A - A**2, etc.
R2 = lambda A: 0.0*A
R3 = lambda A: 0.0*A

# --- Problem (NLBVP) ---
problem = d3.NLBVP([c1, c2, c3, tau1, tau2, tau3], namespace=locals())

# PDEs in each patch: -D*lap(c) + R(c) + tau-lift term = 0
# use lift of tau to match tau polynomial space (see tau method docs)
lift1 = lambda A: d3.Lift(A, basis1, -1)
lift2 = lambda A: d3.Lift(A, basis2, -1)
lift3 = lambda A: d3.Lift(A, basis3, -1)

problem.add_equation(" - {D1} * (d3.Differentiate(c1, r1coord, 2) + (2/r1) * d3.Differentiate(c1, r1coord)) + lift1(tau1) = - (R1(c1)) ".format(D1=D1))
problem.add_equation(" - {D2} * (d3.Differentiate(c2, r2coord, 2) + (2/r2) * d3.Differentiate(c2, r2coord)) + lift2(tau2) = - (R2(c2)) ".format(D2=D2))
problem.add_equation(" - {D3} * (d3.Differentiate(c3, r3coord, 2) + (2/r3) * d3.Differentiate(c3, r3coord)) + lift3(tau3) = - (R3(c3)) ".format(D3=D3))

# --- boundary / interface conditions ---
# center regularity at r=r0 (left of domain1): dr(c1)(r='left') = 0
problem.add_equation("d3.Differentiate(c1, r1coord)(r1='left') = 0")

# outer boundary at r=r3: example Dirichlet c3(r='right') = 0 (change as needed)
problem.add_equation("c3(r3='right') = 0")

# interface at r=r1: continuity
problem.add_equation("c1(r1='right') - c2(r2='left') = 0")

# flux jump across membrane at r=r1:
# D2 * dr(c2) (right side of left) minus D1 * dr(c1) (left side) = k12*(c2 - c1)
problem.add_equation("{D2}*d3.Differentiate(c2, r2coord)(r2='left') - {D1}*d3.Differentiate(c1, r1coord)(r1='right') = {k}*( c2(r2='left') - c1(r1='right') )".format(D1=D1, D2=D2, k=k12))

# interface at r=r2 similarly
problem.add_equation("c2(r2='right') - c3(r3='left') = 0")
problem.add_equation("{D3}*d3.Differentiate(c3, r3coord)(r3='left') - {D2}*d3.Differentiate(c2, r2coord)(r2='right') = {k}*( c3(r3='left') - c2(r2='right') )".format(D2=D2, D3=D3, k=k23))

# --- build solver and Newton iterations (NLBVP) ---
solver = problem.build_solver(ncc_cutoff=1e-12)

# initial guesses (important for Newton)
r1g = dist1.local_grids(basis1)[0]
r2g = dist2.local_grids(basis2)[0]
r3g = dist3.local_grids(basis3)[0]
c1['g'] = 1.0 - (r1g - r0)/(r3 - r0)    # simple ramp guess
c2['g'] = 1.0 - (r2g - r0)/(r3 - r0)
c3['g'] = 1.0 - (r3g - r0)/(r3 - r0)

pert_norm = np.inf
while pert_norm > tolerance:
    solver.newton_iteration()
    # perturbations stored in solver.perturbations; sum norms (as in Lane-Emden example)
    pert_norm = sum(p.allreduce_data_norm('c', 2) for p in solver.perturbations)
    print("perturbation norm:", pert_norm)

# solution is in c1['g'], c2['g'], c3['g']
# %%
import numpy as np
from scipy.sparse import lil_matrix, csr_matrix
from scipy.sparse.linalg import spsolve

# -----------------------
# Problem setup
# -----------------------
a, b = 0.0, 1.0
N = 101                        # number of grid points (uniform grid)
x = np.linspace(a, b, N)
h = x[1] - x[0]

# interface
c = 0.5
ic = np.argmin(np.abs(x - c))   # grid index nearest to c
if ic == 0 or ic == N-1:
    raise ValueError("interface must be strictly internal to grid")

# diffusivities: u has left and right D
Du_left  = 1.0
Du_right = 0.2   # different material to the right
Dv = 0.5        # v diffusivity (continuous)

# partition coefficient for u: u(c^+) = K * u(c^-)
K = 2.0     # set K=1.0 for continuity of u, otherwise partitioning

# flux jumps: here zero because flux continuous
Ju = 0.0
Jv = 0.0

# Dirichlet BCs (example)
u_a, u_b = 0.0, 1.0
v_a, v_b = 0.0, 0.0

# -----------------------
# Reaction (nonlinear) example
# -----------------------
def f(u, v):
    return u - u**3 - v

def g(u, v):
    eps = 0.1
    return eps*(u - 0.5*v)

def fu(u, v): return 1 - 3*u**2
def fv(u, v): return -1.0
def gu(u, v): return 0.1
def gv(u, v): return -0.05

# -----------------------
# Unknown indexing with duplication for u at interface
# -----------------------
# For u, we keep:
#   left side nodes: indices 0 ... ic  (ic is left interface node u_{c^-})
#   right side nodes: indices ic ... N-1 (ic is right interface node u_{c^+})
# This produces N + 1 unknowns for u (duplicate at ic).
Nu = N + 1   # number of u unknowns (duplicate)
Nv = N       # number of v unknowns (continuous)

# Build mapping functions:
# u left indices -> 0..ic
def idx_u_left(i):   # i is grid index 0..ic
    return i
# u right indices -> ic..N (we shift by +ic so that u_right(ic) maps to ic and u_right(ic+1) maps to ic+1 ... up to N)
def idx_u_right(i):  # i is grid index ic..N-1
    # map grid index i to u-index: left side used 0..ic, right uses ic..N-1 -> but we need one extra slot
    # We'll assign u indices 0..ic (left), ic..N (right) => total N+1 entries indexed 0..N
    return i + 0    # this will collide at ic, so we will use convention below:
# To avoid confusion, do mapping explicitly below

# Simpler approach: create a full u_index array of length N+1 mapping each "u-unknown index" to grid position:
# u_unknown indices: 0..N   (N+1 entries)
# represent them as u_unknown_to_grid:
#   for unknown index k in 0..N:
#     k=0..ic   -> grid 0..ic (left)
#     k=ic+1..N -> grid ic..N-1 (right), with k=ic+1 -> grid ic (right), etc.
u_unknown_to_grid = []
# left side unknowns (k=0..ic) map to grid 0..ic
for i in range(0, ic+1):
    u_unknown_to_grid.append(i)
# right side unknowns (k=ic+1..N) map to grid ic..N-1
for i in range(ic, N):
    u_unknown_to_grid.append(i)
# Now len(u_unknown_to_grid) == N+1
# Build reverse mapping: for each grid index and side, where is the u-unknown index?
# We will need to find the u-unknown index for:
#  - left node at grid i (0..ic) -> find smallest k with u_unknown_to_grid[k] == i (the left one)
#  - right node at grid i (ic..N-1) -> find largest k with u_unknown_to_grid[k] == i (the right one)
u_left_idx = np.full(N, -1, dtype=int)   # for grid indices 0..ic -> left u-index
u_right_idx = np.full(N, -1, dtype=int)  # for grid indices ic..N-1 -> right u-index
for k, gidx in enumerate(u_unknown_to_grid):
    # first occurrence is left mapping if multiple occurrences
    if u_left_idx[gidx] == -1:
        u_left_idx[gidx] = k
    # always overwrite right mapping to get the right-side (later occurrence) mapping
    u_right_idx[gidx] = k

# v mapping: continuous 0..N-1, but v unknown indices follow after u unknowns
def idx_v(i):
    return Nu + i   # i = 0..N-1

# u index accessor: specify side 'L' or 'R' for grid index i
def idx_u(i, side='L'):
    if side == 'L':
        return u_left_idx[i]
    else:
        return u_right_idx[i]

# Total unknowns
M = Nu + Nv

# -----------------------
# Residual and Jacobian assembly
# -----------------------
def assemble_res_and_jac(y):
    R = np.zeros(M)
    J = lil_matrix((M, M))

    # unpack u (length Nu) and v (length Nv)
    u = y[0:Nu].copy()
    v = y[Nu:Nu+Nv].copy()

    # Helper to set BC rows for u and v
    # --- u interior on left subdomain: grid i = 0..ic-1 (internal indices)
    for i in range(1, ic):
        ku = idx_u(i, 'L')
        # finite difference using u_left indices
        k_im1 = idx_u(i-1, 'L')
        k_ip1 = idx_u(i+1, 'L')
        # diffusion coefficient on left
        R[ku] = -Du_left * (u[k_im1] - 2*u[ku] + u[k_ip1]) / h**2 + f(u[ku], v[i])
        J[ku, k_im1] = -Du_left / h**2
        J[ku, ku]   =  2*Du_left / h**2 + fu(u[ku], v[i])
        J[ku, k_ip1] = -Du_left / h**2
        J[ku, Nu + i] = fv(u[ku], v[i])

    # --- u interior on right subdomain: grid i = ic+1..N-2
    for i in range(ic+1, N-1):
        ku = idx_u(i, 'R')
        k_im1 = idx_u(i-1, 'R')
        k_ip1 = idx_u(i+1, 'R')
        R[ku] = -Du_right * (u[k_im1] - 2*u[ku] + u[k_ip1]) / h**2 + f(u[ku], v[i])
        J[ku, k_im1] = -Du_right / h**2
        J[ku, ku]   =  2*Du_right / h**2 + fu(u[ku], v[i])
        J[ku, k_ip1] = -Du_right / h**2
        J[ku, Nu + i] = fv(u[ku], v[i])

    # --- v-equations interior for i=1..N-2 (continuous)
    for i in range(1, N-1):
        kv = idx_v(i)
        R[kv] = -Dv * (v[i-1] - 2*v[i] + v[i+1]) / h**2 + g(u[idx_u(i,'L')], v[i])
        J[kv, idx_v(i-1)] = -Dv / h**2
        J[kv, idx_v(i)]   =  2*Dv / h**2 + gv(u[idx_u(i,'L')], v[i])
        J[kv, idx_v(i+1)] = -Dv / h**2
        # coupling to u: use left-side u unknown at grid i (for reaction). 
        # If reaction depends on right-side u on the right subdomain, adjust accordingly.
        J[kv, idx_u(i, 'L')] = gu(u[idx_u(i,'L')], v[i])

    # --- boundary conditions
    # u left boundary: grid 0 -> u_left index at grid 0
    ku0 = idx_u(0, 'L')
    R[ku0] = u[ku0] - u_a
    J[ku0, ku0] = 1.0
    # u right boundary: grid N-1 -> u_right index at grid N-1
    kun = idx_u(N-1, 'R')
    R[kun] = u[kun] - u_b
    J[kun, kun] = 1.0

    # v boundaries
    kv0 = idx_v(0)
    R[kv0] = v[0] - v_a
    J[kv0, kv0] = 1.0
    kvn = idx_v(N-1)
    R[kvn] = v[N-1] - v_b
    J[kvn, kvn] = 1.0

    # --- interface conditions at grid index ic
    # 1) partition: u(c^+) - K * u(c^-) = 0
    kuL = idx_u(ic, 'L')   # left-side u unknown at interface
    kuR = idx_u(ic, 'R')   # right-side u unknown at interface
    R[kuR] = u[kuR] - K * u[kuL]
    J[kuR, kuR] = 1.0
    J[kuR, kuL] = -K

    # 2) flux continuity: Du_left * (u(c^-) - u(c^- - 1))/h  - Du_right * (u(c^+ + 1) - u(c^+))/h = 0
    # left derivative: (u_c_minus - u_{c-1})/h ; right derivative: (u_{c+1} - u_c_plus)/h
    # Use residual row at kuL (or new row). We'll use a separate equation slot — choose kuL to overwrite.
    # Build flux residual R_flux and overwrite kuL row
    kuL = idx_u(ic, 'L')
    u_im1 = idx_u(ic-1, 'L')
    u_ip1 = idx_u(ic+1, 'R')
    # Left derivative: (u[kuL] - u[u_im1]) / h ; Right derivative: (u[u_ip1] - u[kuR]) / h
    R_flux = Du_left * (u[kuL] - u[u_im1]) / h - Du_right * (u[u_ip1] - u[kuR]) / h - Ju
    R[kuL] = R_flux
    # Jacobian entries for flux eqn
    J[kuL, kuL] = Du_left / h
    J[kuL, u_im1] = -Du_left / h
    J[kuL, u_ip1] = -Du_right / h
    J[kuL, kuR]   = Du_right / h

    # 3) v at interface: we kept v continuous, so v-equation at i=ic must be the usual FD (we already did it above)
    # But ensure its Jacobian used correct u variable (here use left u mapping for reaction)
    # (Already set above in v interior loop for i=ic)

    return R, J

# -----------------------
# Newton solver
# -----------------------
def newton(y0, tol=1e-8, maxit=40):
    y = y0.copy()
    for it in range(maxit):
        R, Jl = assemble_res_and_jac(y)
        normR = np.linalg.norm(R, np.inf)
        print(f"iter {it}, ||R||_inf = {normR:.3e}")
        if normR < tol:
            return y, True
        J = csr_matrix(Jl)
        try:
            dy = spsolve(J, -R)
        except Exception as e:
            print("Linear solve failed:", e)
            return y, False
        # simple damping
        alpha = 1.0
        y_new = y + alpha*dy
        for ls in range(10):
            Rn, _ = assemble_res_and_jac(y_new)
            if np.linalg.norm(Rn, np.inf) < normR:
                break
            alpha *= 0.5
            y_new = y + alpha*dy
        y = y_new
    return y, False

# -----------------------
# initial guess
# -----------------------
# u initial: linear left and right with jump ratio K at interface
u_init = np.zeros(Nu)
# left u unknowns k=0..ic map to grid 0..ic : linear from u_a to some mid value
for k in range(0, ic+1):
    g = u_unknown_to_grid[k]
    # simple linear interpolation from a to c using grid position
    u_init[k] = u_a + (u_b - u_a)*(x[g] - a)/(b - a) * 0.8

# right unknowns k=ic+1..N map to grid ic..N-1
for k in range(ic+1, Nu):
    g = u_unknown_to_grid[k]
    u_init[k] = u_a + (u_b - u_a)*(x[g] - a)/(b - a) * 0.8

# enforce partition approximately: set kuR = K * kuL
u_init[idx_u(ic,'R')] = K * u_init[idx_u(ic,'L')]

v_init = np.linspace(v_a, v_b, N)

y0 = np.concatenate([u_init, v_init])

# solve
y_sol, ok = newton(y0)
if not ok:
    print("Newton failed to converge.")
else:
    print("Converged.")
    u_sol = y_sol[0:Nu]
    v_sol = y_sol[Nu:Nu+Nv]
    # map u unknowns back to left/right grid values for plotting if desired
    u_left_vals  = np.array([u_sol[idx_u(i,'L')] for i in range(0, ic+1)])
    u_right_vals = np.array([u_sol[idx_u(i,'R')] for i in range(ic, N)])
    # construct a piecewise plot vector
    x_u = np.concatenate([x[0:ic+1], x[ic:N]])
    u_plot = np.concatenate([u_left_vals, u_right_vals])

    try:
        import matplotlib.pyplot as plt
        plt.plot(x_u, u_plot, marker='o', label='u (dup nodes)')
        plt.plot(x, v_sol, marker='x', label='v (continuous)')
        plt.axvline(x[ic], color='k', linestyle='--', label='interface')
        plt.legend()
        plt.show()
    except Exception:
        pass

# %%
