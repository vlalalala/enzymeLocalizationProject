import numpy as np
from scipy.integrate import solve_ivp
from dataclasses import dataclass

class SystemParams:
    """
    Parameters defining the reaction-diffusion system.

    Attributes
    ----------
    radii : array (N,)
        Outer radii R_1 < ... < R_N. Innermost compartment starts at r=0.
    D : array (M,)
        Diffusion coefficients per species, constant across compartments.
    K : list of N arrays, each (M, M)
        K[n] is the reaction rate matrix in compartment n.
        For a linear chain X_1->X_2->...->X_M:
            K[n][i,i]   = -k_i^(n)   (diagonal: loss from species i)
            K[n][i+1,i] = +k_i^(n)   (subdiagonal: gain to species i+1)
    P : array (N-1, M)
        P[n, i] = permeability of species i at interface n.
    P_out : array (M,)
        Permeability of each species at the outer boundary.
    q_inf : array (M,)
        External concentrations at the outer boundary.
    """

    def __init__(self, radii, D, K, P, P_out, q_inf):
        self.radii = radii #np.asarray(self.radii, dtype=float) # (N,)
        self.D     = D #np.asarray(self.D,     dtype=float) # (M,)
        self.K     = K #[np.asarray(Kn, dtype=float) for Kn in self.K] # N arrays of shape (M, M)
        self.P     = P #np.asarray(self.P,     dtype=float) # (N-1, M)
        self.P_out = P_out #np.asarray(self.P_out, dtype=float) # (M,)
        self.q_inf = q_inf #np.asarray(self.q_inf, dtype=float) # (M,)
        self.N = len(self.radii) # number of compartments
        self.M = len(self.D) # number of species
        if self.N > 1:
            assert self.P.shape == (self.N - 1, self.M), \
                f"P must have shape (N-1, M) = ({self.N-1}, {self.M})"
        assert len(self.K) == self.N
        for n, Kn in enumerate(self.K):
            assert Kn.shape == (self.M, self.M), \
                f"K[{n}] must have shape (M, M) = ({self.M}, {self.M})"

def nondimensionalize(params):
    """
    Rescale to dimensionless units for numerical stability.

    Scales:
      Length : L  = R_N
      Diff.  : D0 = D[0]
      Conc.  : C  = max(|q_inf|), or 1 if all zero

    Dimensionless groups:
      r_tilde = r / L
      q_tilde = q / C
      J_tilde = J * L / (D0 * C)     [from J = -D dq/dr]
      K_tilde = K * L^2 / D0
      P_tilde = P * L / D0
      D_tilde = D / D0
    """
    L  = params.radii[-1]
    D0 = params.D[0]
    C  = np.max(np.abs(params.q_inf))
    if C == 0:
        C = 1.0

    nd = SystemParams(
        radii = params.radii / L,
        D     = params.D     / D0,
        K     = [Kn * (L**2 / D0) for Kn in params.K],
        P     = params.P     * (L / D0) if params.N > 1 else params.P,
        P_out = params.P_out * (L / D0),
        q_inf = params.q_inf / C,
    )
    return nd, {'L': L, 'D0': D0, 'C': C}

def make_A(r, Kn, D_inv):
    """
    Build the 2M x 2M system matrix at radius r for compartment n.

    A(r) = [ 0       -D^{-1} ]
           [ K^(n)   -2/r I  ]

    State vector z = (q, J) in R^{2M}.
    """
    M = len(D_inv) # number of species
    A = np.zeros((2*M, 2*M))
    A[:M, M:] = -np.diag(D_inv)   # top-right: -D^{-1}
    A[M:, :M] = Kn                 # bottom-left: K^(n)
    A[M:, M:] = -(2.0/r) * np.eye(M)  # bottom-right: -2/r I
    return A

def compute_transfer_matrix(r_start, r_end, Kn, D_inv, rtol=1e-10, atol=1e-12):
    """
    Compute T_n in R^{2M x 2M} by integrating the matrix ODE

        dT/dr = A^(n)(r) T,    T(r_start) = I_{2M}

    from r_start to r_end. Each column of T is integrated as a
    separate ODE solution with initial condition e_j.

    Parameters
    ----------
    r_start : float  (= R_{n-1}, dimensionless)
    r_end   : float  (= R_n,     dimensionless)
    Kn      : (M, M) reaction matrix for compartment n (dimensionless)
    D_inv   : (M,)   1/D_i (dimensionless)
    rtol, atol : ODE solver tolerances

    Returns
    -------
    T : (2M, 2M) transfer matrix
    """
    M2 = 2 * len(D_inv)

    def rhs(r, T_flat):
        T = T_flat.reshape(M2, M2)
        A = make_A(r, Kn, D_inv)
        return (A @ T).ravel()

    T0 = np.eye(M2).ravel()
    sol = solve_ivp(
        rhs,
        [r_start, r_end],
        T0,
        #method="Radau",
        method='DOP853',
        rtol=rtol,
        atol=atol,
        dense_output=False,
    )
    if not sol.success:
        raise RuntimeError(
            f"ODE integration failed for compartment [{r_start:.4f}, {r_end:.4f}]: "
            + sol.message
        )
    return sol.y[:, -1].reshape(M2, M2)


def interface_matrix(P_n, M):
    """
    Build the 2M x 2M interface matrix at an internal interface.

    Permeability condition: J^- = P_n (q^- - q^+)
    =>  q^+ = q^- - P_n^{-1} J^-
        J^+ = J^-

    So: I_n = [ I   -P_n^{-1} ]
              [ 0       I     ]

    Parameters
    ----------
    P_n : (M,) permeabilities at this interface (dimensionless)
    M   : number of species

    Returns
    -------
    In : (2M, 2M) interface matrix
    """
    In = np.eye(2 * M)
    In[:M, M:] = -np.diag(1.0 / P_n)
    return In


# ─────────────────────────────────────────────────────────────
# Global propagator G
# ─────────────────────────────────────────────────────────────

def compute_global_propagator(nd):
    """
    Compute G = T_N I_{N-1} T_{N-1} ... I_1 T_1 in R^{2M x 2M}.

    Applies right-to-left: innermost compartment first.
    """
    M  = nd.M
    N  = nd.N
    D_inv = 1.0 / nd.D

    # start with innermost compartment
    # r_start for compartment 1 is r=0, but A(r) is singular at r=0.
    # Use a small epsilon instead; the regularity BC handles r=0 separately.
    r_inner = nd.radii[0]
    eps = r_inner * 1e-6
    G = compute_transfer_matrix(eps, r_inner, nd.K[0], D_inv)

    for n in range(1, N):
        r_left  = nd.radii[n-1]
        r_right = nd.radii[n]

        # interface matrix at r = R_{n-1} (between compartments n-1 and n)
        In = interface_matrix(nd.P[n-1], M)
        G  = compute_transfer_matrix(r_left, r_right, nd.K[n], D_inv) @ In @ G

    return G


# ─────────────────────────────────────────────────────────────
# Solve the system
# ─────────────────────────────────────────────────────────────

def solve(params: SystemParams):
    """
    Solve for the steady-state concentration profiles.

    Steps:
      1. Nondimensionalize
      2. Build global propagator G
      3. Apply inner BC (B_0) and outer Robin BC (C_out)
      4. Solve m x m system (C_out G B_0) a = b for a in R^m
      5. Return solution object for evaluation

    Returns
    -------
    sol : dict with fields needed by evaluate_solution()
    """
    nd, scales = nondimensionalize(params)
    M = nd.M

    # Inner BC: z(0^+) = B_0 a,  B_0 = [I; 0] in R^{2M x M}
    eps = nd.radii[0] * 1e-6
    B0 = np.zeros((2*M, M))
    B0[:M, :] = np.eye(M)
    B0[M:, :] = (eps / 3.0) * nd.K[0] # Matches the regular J profile at r=eps ########################################
    # Global propagator
    G = compute_global_propagator(nd)

    # Outer Robin BC: J + P_out q = P_out q_inf
    # => C_out z = b,  C_out = [P_out | I] in R^{M x 2M}, b = P_out q_inf
    P_out_mat = np.diag(nd.P_out)
    C_out = np.zeros((M, 2*M))
    C_out[:, :M] = P_out_mat   # coefficient of q
    C_out[:, M:] = -np.eye(M)   # coefficient of J ############################################################

    b = nd.P_out * nd.q_inf

    # Reduced system: (C_out G B_0) a = b
    A_sys = C_out @ G @ B0   # M x M
    a     = np.linalg.solve(A_sys, b)

    return {
        'a':      a,
        'G':      G,
        'B0':     B0,
        'nd':     nd,
        'scales': scales,
        'params': params,
    }


# ─────────────────────────────────────────────────────────────
# Evaluate solution on a radial grid
# ─────────────────────────────────────────────────────────────

def evaluate_solution(sol, r_points):
    """
    Evaluate concentration profiles at dimensional radii r_points.

    Propagates z from r=0 outward, tracking which compartment we are in,
    and integrating the ODE up to each requested r.

    Returns
    -------
    q_out : array (M, len(r_points))
        q_out[i, k] = concentration of species i at r_points[k]
    """
    params  = sol['params']
    nd      = sol['nd']
    scales  = sol['scales']
    a       = sol['a']
    M       = params.M
    N       = params.N
    L       = scales['L']
    C       = scales['C']
    D_inv   = 1.0 / nd.D

    # Ensure points are sorted to allow single-pass integration
    r_points = np.sort(np.asarray(r_points, dtype=float))
    q_out    = np.zeros((M, len(r_points)))

    z = sol['B0'] @ a   # Start at r = eps
    eps = nd.radii[0] * 1e-6
    r_prev = eps
    n = 0 # current compartment tracker

    for k, r_dim in enumerate(r_points):
        r_nd = r_dim / L
        
        # Move through whole compartments if r_nd is further out
        while n < N and r_nd > nd.radii[n]:
            if nd.radii[n] > r_prev:
                T = compute_transfer_matrix(r_prev, nd.radii[n], nd.K[n], D_inv)
                z = T @ z
            
            # Apply internal interface map
            if n < N - 1:
                In = interface_matrix(nd.P[n], M)
                z = In @ z
                
            r_prev = nd.radii[n]
            n += 1

        # Integrate the remainder of the current compartment up to r_nd
        if r_nd > r_prev:
            T = compute_transfer_matrix(r_prev, r_nd, nd.K[n], D_inv)
            z = T @ z
            r_prev = r_nd

        q_out[:, k] = C * z[:M]

    return q_out


# ─────────────────────────────────────────────────────────────
# Helper: build K matrix for a linear chain
# ─────────────────────────────────────────────────────────────

def build_K_chain(k_rates):
    """
    Build the M x M reaction matrix K for a linear chain
    X_1 -> X_2 -> ... -> X_M with rates k_rates = [k_1, ..., k_{M-1}].

    K[i,i]   = -k_i     (diagonal: loss from species i, for i < M-1)
    K[i+1,i] = +k_i     (subdiagonal: gain to species i+1)
    Last species has no sink: K[M-1, M-1] = 0
    """
    k = np.asarray(k_rates, dtype=float)
    M = len(k) + 1
    K = np.zeros((M, M))
    for i in range(M - 1):
        K[i,   i] = -k[i]
        K[i+1, i] = +k[i]
    return K


# ─────────────────────────────────────────────────────────────
# Example
# ─────────────────────────────────────────────────────────────

if __name__ == '__main__':
    import matplotlib.pyplot as plt

    r_plot = np.linspace(0.01, 3.0, 500)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))

    # ── 2-species, 3-compartment ──
    N, M = 3, 2
    params2 = SystemParams(
        radii = np.array([1.0, 2.0, 3.0]),
        D     = np.array([1.0, 1.0]),
        K     = [build_K_chain([k]) for k in [0.5, 1.0, 2.0]],
        P     = np.array([[1.0, 1.0],
                          [1.0, 1.0]]),
        P_out = np.array([1.0, 1.0]),
        q_inf = np.array([1.0, 0.0]),
    )
    sol2 = solve(params2)
    q2   = evaluate_solution(sol2, r_plot)

    ax = axes[0]
    ax.plot(r_plot, q2[0], label='$X_1$')
    ax.plot(r_plot, q2[1], label='$X_2$')
    for R in params2.radii[:-1]:
        ax.axvline(R, color='gray', linestyle='--', lw=0.8)
    ax.set_xlabel('r'); ax.set_ylabel('Concentration')
    ax.set_title('2-species chain, 3 compartments')
    ax.legend()

    # ── 3-species, 3-compartment ──
    params3 = SystemParams(
        radii = np.array([1.0, 2.0, 3.0]),
        D     = np.array([1.0, 1.0, 1.0]),
        K     = [build_K_chain([k1, k2]) for k1, k2 in [(0.5,0.3),(1.0,0.6),(2.0,1.2)]],
        P     = np.array([[1.0, 1.0, 1.0],
                          [1.0, 1.0, 1.0]]),
        P_out = np.array([1.0, 1.0, 1.0]),
        q_inf = np.array([1.0, 0.0, 0.0]),
    )
    sol3 = solve(params3)
    q3   = evaluate_solution(sol3, r_plot)

    ax = axes[1]
    ax.plot(r_plot, q3[0], label='$X_1$')
    ax.plot(r_plot, q3[1], label='$X_2$')
    ax.plot(r_plot, q3[2], label='$X_3$')
    for R in params3.radii[:-1]:
        ax.axvline(R, color='gray', linestyle='--', lw=0.8)
    ax.set_xlabel('r'); ax.set_ylabel('Concentration')
    ax.set_title('3-species chain, 3 compartments')
    ax.legend()

    plt.tight_layout()
    plt.savefig('/mnt/user-data/outputs/reaction_diffusion_profiles.png', dpi=150)
    print('Done.')




