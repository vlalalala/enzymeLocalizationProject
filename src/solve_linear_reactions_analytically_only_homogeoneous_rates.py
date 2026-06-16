import numpy as np

def calculate_analytical_solution_1_spontaneous_reaction_1_region(k, D, R, p, ext):
    """Returns a lambda function
    D is diffusion coefficient for reactant, p is permeability constant of reactant,
    ext is external concentration of reactant
    """
    s_lambda = np.sqrt(k / D)
    A = - ((p / D * ext * R**2)
        / (
            np.exp(s_lambda * R)*(s_lambda * R - 1 + (p * R)/D)
            + np.exp(-s_lambda * R) * (s_lambda * R + 1 - (p * R)/D) 
            )
    )
    c = lambda r: 1/r * A *(np.exp(-s_lambda*r)-np.exp(s_lambda*r))
    return c

def calculate_analytical_solution_1_spontaneous_reaction_2_regions(k, D, R, p, ext, r_inner):
    X_lambda = np.sqrt(k / D)
    s = np.sinh(X_lambda * r_inner)
    c = np.cosh(X_lambda * r_inner)
    beta = p / D
    alpha = p * R / D
    
    rho = np.exp(-2 * X_lambda * r_inner) * (
        D * (r_inner**2 * X_lambda**2 * c + r_inner * X_lambda * (c - s) - s)
        + p * r_inner**2 * X_lambda * (s + c)
    ) / (
        D * (r_inner**2 * X_lambda**2 * c - r_inner * X_lambda * (s + c) + s)
        + p * r_inner**2 * X_lambda * (s - c)
    )

    A = (beta * ext * R**2)/(
        np.exp(-X_lambda * R)*(alpha - X_lambda * R - 1) + rho * np.exp(X_lambda * R)*(alpha + X_lambda*R -1)
    )
    
    B = rho * A

    S = (
        beta * r_inner * (A * np.exp(-X_lambda * r_inner) + B * np.exp(X_lambda * r_inner))
    ) / (
        X_lambda * r_inner * np.cosh(X_lambda * r_inner)
        - np.sinh(X_lambda * r_inner)
        + beta * r_inner * np.sinh(X_lambda * r_inner)
    )

    
    c_1 = lambda r : S * np.sinh(X_lambda * r) / r
    c_2 = lambda r : (A * np.exp(-X_lambda * r) + B * np.exp(X_lambda * r))/r
    c = lambda r: np.where(
        r < r_inner,
        c_1(r),
        c_2(r)
    )
    return c


def add_theory_curve_to_ax(
        fig,
        ax,
        reaction_network,
        num_regions,
        membrane_radii,
    ):
    membrane_radii_with_0 = [0] + membrane_radii

    # Case of one spontaneous reaction with no inner boundaries
    if (len(reaction_network.spontaneous_reactions) == 1
        and len(reaction_network.enzymatic_reactions) == 0
        and num_regions == 1
    ):  
        s_reaction = reaction_network.spontaneous_reactions[0]
        s = s_reaction.start_species
        external_radius = membrane_radii[-1]
        c = calculate_analytical_solution_1_spontaneous_reaction_1_region(
            k = s_reaction.k,
            D = s.diffusion_constant,
            R = external_radius,
            p = s.permeability_constant,
            ext = s.external_concentration
        )
        """
        print("Plotting analytical solution")
        s_reaction = reaction_network.spontaneous_reactions[0]
        s = s_reaction.start_species
        s_lambda = np.sqrt(s_reaction.k / sion_constant)
        external_radius = membrane_radii[-1]
        A = - ((s.permeability_constant / s.diffusion_constant * s.external_concentration * external_radius**2)
            / (
                np.exp(s_lambda * external_radius)*(s_lambda * external_radius - 1 + (s.permeability_constant * external_radius)/s.diffusion_constant)
               + np.exp(-s_lambda * external_radius) * (s_lambda * external_radius + 1 - (s.permeability_constant * external_radius)/s.diffusion_constant) 
              )
        )
        c = lambda r: 1/r * A *(np.exp(-s_lambda*r)-np.exp(s_lambda*r))
        """
        r_to_plot = np.linspace(external_radius*0.01, external_radius, num = 100)
        ax.plot([r/external_radius for r in r_to_plot], [c(r) for r in r_to_plot], linestyle= "--",
                label = f"analytical solution for {s.name}", 
                linewidth = 1,
                color = "k",
                alpha = 0.5,
                zorder = 100
        )
        ax.legend()
    
    # Case of one spontaneous reaction with one inner boundary
    elif (len(reaction_network.spontaneous_reactions) == 1
        and len(reaction_network.enzymatic_reactions) == 0
        and num_regions == 2
    ):  
        reaction = reaction_network.spontaneous_reactions[0]
        X = reaction.start_species
        external_radius = membrane_radii[-1]
        c = calculate_analytical_solution_1_spontaneous_reaction_2_regions(
            k = reaction.k,
            D = X.diffusion_constant,
            R = external_radius,
            p = X.permeability_constant,
            ext = X.external_concentration,
            r_inner = membrane_radii[0]
        )
        print("Plotting analytical solution")
        """
        X_lambda = np.sqrt(reaction.k / X.diffusion_constant)
        r_inner = membrane_radii[0]
        external_radius = membrane_radii[-1]
        s = np.sinh(X_lambda * r_inner)
        c = np.cosh(X_lambda * r_inner)
        beta = X.permeability_constant / X.diffusion_constant
        alpha = X.permeability_constant * external_radius / X.diffusion_constant
        
        rho = np.exp(-2 * X_lambda * r_inner) * (
            X.diffusion_constant * (r_inner**2 * X_lambda**2 * c + r_inner * X_lambda * (c - s) - s)
            + X.permeability_constant * r_inner**2 * X_lambda * (s + c)
        ) / (
            X.diffusion_constant * (r_inner**2 * X_lambda**2 * c - r_inner * X_lambda * (s + c) + s)
            + X.permeability_constant * r_inner**2 * X_lambda * (s - c)
        )

        A = (beta * X.external_concentration * external_radius**2)/(
            np.exp(-X_lambda * external_radius)*(alpha - X_lambda * external_radius - 1) + rho * np.exp(X_lambda * external_radius)*(alpha + X_lambda*external_radius -1)
        )
        
        B = rho * A

        S = (
            beta * r_inner * (A * np.exp(-X_lambda * r_inner) + B * np.exp(X_lambda * r_inner))
        ) / (
            X_lambda * r_inner * np.cosh(X_lambda * r_inner)
            - np.sinh(X_lambda * r_inner)
            + beta * r_inner * np.sinh(X_lambda * r_inner)
        )

        
        c_1 = lambda r : S * np.sinh(X_lambda * r) / r
        c_2 = lambda r : (A * np.exp(-X_lambda * r) + B * np.exp(X_lambda * r))/r
        c = [c_1, c_2]
        """
        for region_idx in [0,1]:
            max_radius = membrane_radii_with_0[region_idx+1]
            min_radius = membrane_radii_with_0[region_idx]
            region_radii = np.linspace(min_radius, max_radius)
            # region_radii[0] skipped to avoid division by 0
            r_to_plot = np.linspace(region_radii[0]+external_radius*0.01, region_radii[-1], num = 100)
            if region_idx == 0:
                label = f"theory for {X.name}"
            else:
                label = None
            ax.plot([r/external_radius for r in r_to_plot],
                [c(r) for r in r_to_plot],
                #[c[region_idx](r) for r in r_to_plot],
                linestyle= "--",
                label = label, zorder = 100, 
                linewidth = 1,
                color = "k",
                alpha = 0.5
            )
        ax.legend()
    
    # Case of one spontaneous reaction with two inner boundaries
    elif (len(reaction_network.spontaneous_reactions) == 1
        and len(reaction_network.enzymatic_reactions) == 0
        and num_regions == 3
    ):  
        print("Plotting analytical solution")

        """This below was written by ChatGPT."""
        reaction = reaction_network.spontaneous_reactions[0]
        X = reaction.start_species

        lam = np.sqrt(reaction.k / X.diffusion_constant)  # lambda
        D = X.diffusion_constant
        p = X.permeability_constant
        beta = p / D

        # Two inner membranes
        R1 = membrane_radii[0]
        R2 = membrane_radii[1]
        R  = membrane_radii[2]

        c_ext = X.external_concentration

        # --- Helpers for evaluating c and c' for the basis solutions ---
        # For region 1: c1(r) = S1*sinh(lam r)/r
        def c1_val(r, S1):
            return S1 * np.sinh(lam * r) / r

        def c1_der(r, S1):
            # d/dr [ S1*sinh(lam r)/r ] = S1*(lam*cosh(lam r)/r - sinh(lam r)/r^2)
            return S1 * (lam * np.cosh(lam * r) / r - np.sinh(lam * r) / (r**2))

        # For region j>=2: cj(r) = (A*e^{-lam r} + B*e^{lam r})/r
        def cAB_val(r, A, B):
            return (A * np.exp(-lam * r) + B * np.exp(lam * r)) / r

        def cAB_der(r, A, B):
            # derivative of (f(r)/r) with f=A e^{-lam r}+B e^{lam r}
            f  = A * np.exp(-lam * r) + B * np.exp(lam * r)
            fp = -lam * A * np.exp(-lam * r) + lam * B * np.exp(lam * r)
            return fp / r - f / (r**2)

        # --- Build the 5x5 linear system M x = b for x=[S1, A2, B2, A3, B3] ---
        # Conditions:
        # (1) c1'(R1) = c2'(R1)
        # (2) c1'(R1) = beta*(c2(R1) - c1(R1))
        # (3) c2'(R2) = c3'(R2)
        # (4) c2'(R2) = beta*(c3(R2) - c2(R2))
        # (5) c3'(R)  = beta*(c_ext - c3(R))

        M = np.zeros((5, 5), dtype=float)
        b = np.zeros(5, dtype=float)

        # Row 0: c1'(R1) - c2'(R1) = 0
        # coefficients for S1, A2, B2, A3, B3
        M[0, 0] = c1_der(R1, 1.0)
        # -c2'(R1) contributions:
        # c2'(R1) is linear in A2,B2 so put negatives
        # We'll get coeffs by evaluating derivative with A=1,B=0 and A=0,B=1
        M[0, 1] = -cAB_der(R1, 1.0, 0.0)
        M[0, 2] = -cAB_der(R1, 0.0, 1.0)
        # A3,B3 not in this equation
        b[0] = 0.0

        # Row 1: c1'(R1) - beta*(c2(R1) - c1(R1)) = 0
        # => c1'(R1) + beta*c1(R1) - beta*c2(R1) = 0
        M[1, 0] = c1_der(R1, 1.0) + beta * c1_val(R1, 1.0)
        M[1, 1] = -beta * cAB_val(R1, 1.0, 0.0)
        M[1, 2] = -beta * cAB_val(R1, 0.0, 1.0)
        b[1] = 0.0

        # Row 2: c2'(R2) - c3'(R2) = 0
        M[2, 1] = cAB_der(R2, 1.0, 0.0)
        M[2, 2] = cAB_der(R2, 0.0, 1.0)
        M[2, 3] = -cAB_der(R2, 1.0, 0.0)
        M[2, 4] = -cAB_der(R2, 0.0, 1.0)
        b[2] = 0.0

        # Row 3: c2'(R2) - beta*(c3(R2) - c2(R2)) = 0
        # => c2'(R2) + beta*c2(R2) - beta*c3(R2) = 0
        M[3, 1] = cAB_der(R2, 1.0, 0.0) + beta * cAB_val(R2, 1.0, 0.0)
        M[3, 2] = cAB_der(R2, 0.0, 1.0) + beta * cAB_val(R2, 0.0, 1.0)
        M[3, 3] = -beta * cAB_val(R2, 1.0, 0.0)
        M[3, 4] = -beta * cAB_val(R2, 0.0, 1.0)
        b[3] = 0.0

        # Row 4: c3'(R) = beta*(c_ext - c3(R))
        # => c3'(R) + beta*c3(R) = beta*c_ext
        M[4, 3] = cAB_der(R, 1.0, 0.0) + beta * cAB_val(R, 1.0, 0.0)
        M[4, 4] = cAB_der(R, 0.0, 1.0) + beta * cAB_val(R, 0.0, 1.0)
        b[4] = beta * c_ext

        # Solve
        S1, A2, B2, A3, B3 = np.linalg.solve(M, b)

        # --- Define concentration functions for each region (compatible with your plotting pattern) ---
        def c_1(r):
            r = np.asarray(r, dtype=float)
            out = np.empty_like(r)
            # safe near 0: limit S1*sinh(lam r)/r -> S1*lam
            small = np.isclose(r, 0.0)
            out[small] = S1 * lam
            out[~small] = S1 * np.sinh(lam * r[~small]) / r[~small]
            return out

        def c_2(r):
            r = np.asarray(r, dtype=float)
            return (A2 * np.exp(-lam * r) + B2 * np.exp(lam * r)) / r

        def c_3(r):
            r = np.asarray(r, dtype=float)
            return (A3 * np.exp(-lam * r) + B3 * np.exp(lam * r)) / r

        c = [c_1, c_2, c_3]

        # --- Plot using your mesh_points_in_regions structure (now 3 regions) ---

        for region_idx in [0, 1, 2]:
            max_radius = membrane_radii_with_0[region_idx+1]
            min_radius = membrane_radii_with_0[region_idx]
            region_radii = np.linspace(min_radius, max_radius)
            # skip first point in region 0 to avoid r=0 if present
            start_idx = 1 if (region_idx == 0 and np.isclose(region_radii[0], 0.0)) else 0

            r_to_plot = np.linspace(region_radii[start_idx], region_radii[-1], num=100)

            label = f"theory for {X.name}" if region_idx == 0 else None
            ax.plot([rr / R for rr in r_to_plot],
                    [float(c[region_idx](rr)) for rr in r_to_plot],
                    linestyle="--",
                    label=label,
                    zorder=100,
                    linewidth=1,
                    color = "k",
                    alpha = 0.5
            )
    
    
    
    return fig, ax

       

