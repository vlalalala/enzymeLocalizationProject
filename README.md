# Numerical solver of chemical reactions within concentric spherical semi-permeable membranes

This code is for calculating the steady state concentrations of substances placed in a spherically symmetrical system of semi-permeable membranes, where the concentration of  each of substances outside the exterior membrane is kept constant. Spontaneous and enzymatic reactions can be defined (enzymes can be placed in the regions between semi-permeable membranes) and an arbitrary number of semi-permeable membranes can be used. The boundary problem is solved numerically through the Newton method.
![alt text](examples/simple_decay_with_two_inner_boundaries/combined_000001/newton_iterations.gif)

## Summary

This work continues previous work by Hinzpeter et al
- Optimal Compartmentalization Strategies for Metabolic Microcompartments, by Hinzpeter et al. (Biophysical Journal, 2017)

The governing equations for each species $q$ are reaction-diffusion equations of the form

$$ \frac{\partial q}{\partial t} = D_q \nabla^2 q + R_q $$
with $D_q$ the diffusion constant and $R_q$ the reaction term resulting from interactions between the different species,
and with boundary conditions given by
$$ J = p \cdot (q(r^-) - q(r^+)) $$
where $q(r^-)$ and $q(r^+)$ are the concentrations of $q$ at either side of the boundary. For each species, $q(R^+)$ is constant.

The steady state distribution of $q$ is computed numerically.

The interval $[0, R]$ is discretized by setting equally spaced mesh points. The position of each inner membrane is set to the closest mesh point position. For $N$ inner membranes defined, there are $N+1$ regions. Each region is defined by the mesh points within its bounds, including those at the bounding membranes. (Therefore, at the mesh positions where a inner boundary is at, there are in actuality 2 mesh points.)

The solution is assumed to have converged well enough if the net flux computed (reaction flux - flux through boundary) is close to zero for each species (i.e. mass conservation for each species).

## Theory
In spherical symmetry (assuming no angular dependence), the steady-state reaction-diffusion equation becomes

$$ \frac{1}{r^2} \frac{d}{dr} \left( r^2 \frac{du}{dr} \right) + \frac{1}{D} R(u) = 0 ,  $$

which simplifies to

$$\frac{d^2u}{dr^2} + \frac{2}{r} \frac{du}{dr} + \frac{1}{D} R(u) = 0$$

(Note the singularity at r=0).

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

such that

$$ \left(\frac{dy_0}{dr}(0) = \right)\ \  y_1(0) = 0 $$



and at $r=R$ we have

$$ \frac{du}{dr}(R) = [u_\mathrm{ext} − u(R)] \cdot \frac{p}{D}$$

such that
$$
y_1(R) = (u_{ext} - y_0(R)) \cdot \frac{p}{D}
$$

with $p$ the permeability to the outer membrane for substance $u$.

**Generally**:

For

- $m$ species
- $n$ compartments i.e. $n-1$ inner membranes at $r_i$ with $i \in \{1,2,...,n-1\} $
- $m \cdot [2 + (n-1)\cdot 2] = 2mn$ boundary conditions (for each species: one BC at r=0, one BC at r=R, one BC at each end of each inner membrane)


Define $2\cdot m\cdot n$ variables, with $k\in \{0, 1, ..., n-1\}$

$$v_k^{(m)} = u_k^{(m)}$$
$$ w_k^{(m)} = \frac{\mathrm{d}u_k^{(m)}}{\mathrm{d}r} $$

with 

$$ \frac{\mathrm{d}v_k^{(m)}}{\mathrm{d}r} = w_k^{(m)}$$
$$ \frac{\mathrm{d}w_k^{(m)}}{\mathrm{d}r} = -\frac{2}{r} w_k^{(m)}-\frac{1}{D}R(\vec{v}_k)$$

Boundary conditions are given by

$$ w^{(m)}(0) = 0 $$

(Neumann no-flux BC)

$$ w^{(m)}(R) = \left(u_\mathrm{ext} - v^{(m)}(R)\right) \cdot \frac{p^{(m)}}{D} $$

(Robin BC)

$$w_i^{(m)}(r_i^+) =  \frac{p^{(m)}}{D} \cdot \left( v_{i}^{(m)}(r_i^+) - v_{i-1}^{(m)}(r_i^-)\right)$$

$$ w_{i-1}^{(m)}(r_i^-) =  \frac{p^{(m)}}{D} \cdot \left( v_{i}^{(m)}(r_i^+) - v_{i-1}^{(m)}(r_i^-)\right)$$

The concentration is not continuous at the membranes, but the flux is continuous!

We used

$$f'(x_i) = \frac{1}{2h}(f_{i+1}-f_{i-1}) + O(h^2)$$

$$f''(x_i) \approx \frac{1}{h^2}(f_{i+1}-2f_i + 
f_{i-1})$$

for the discretization.

**Given multiple regions**:

In each region $k$ we have 

$$ \frac{1}{r^2} \frac{d}{dr} \left( r^2 \frac{du_k}{dr} \right) + \frac{1}{D_k} R(u) = 0 ,  $$

$$ D_k \, \frac{\mathrm{d}u_k}{\mathrm{d}r} \big|_{r = R_k^-}
= D_{k+1} \, \frac{\mathrm{d}u_{k+1}}{\mathrm{d}r} \big|_{r = R_k^+}
= p_k \left( u_{k+1}(R_k^+) - u_k(R_k^-) \right)$$

**Fluxes**

**Finding 
To find the total number of moles $n$ within the volume with radius R, calculate
$$
n = \int_0^R C(r) \cdot 4\pi r^2 \mathrm{d}r
$$


Reminder that to convert $C(r)$ from $M$ (moles per litre) to moles per $m^{3}$: $M = \frac{moles}{m^3}\cdot 10^{-3}$.

So 
$$
n = \int_0^R C(r) \cdot 10^{3}\cdot 4\pi r^2 \mathrm{d}r
$$
with $n$ in moles and $C(r)$ in $M$.

## Typical values

Diffusion coefficient of l-Trp taken as $6.6\cdot 10^{-6} cm^{2} s^{-1} = 6.6\cdot 10^{-10} m^{2} s^{-1}$ (estimated from https://pmc.ncbi.nlm.nih.gov/articles/PMC16526/).
Assuming same diffusion coefficient for the other values.

Following https://pubs.acs.org/doi/pdf/10.1021/acsomega.3c08233?ref=article_openPDF (quantification of kinetics of VioA) we assume $k_\mathrm{M} = 125 \mu M = 1.25\cdot 10^{-4}M$ and $k_\mathrm{cat} = 0.75 s^{-1}$ for all enzymes.

For the external concentration of l-Trp we're assuming $s_0 = 25 \mu M = 25\cdot 10^{-6}M$ and the permeability $p = 90 \mu M s^{-1} = 90\cdot 10^{-6}M s^{-1}$.
(as in https://www.cell.com/biophysj/fulltext/S0006-3495(16)34263-1?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0006349516342631%3Fshowall%3Dtrue)

We use a radius of $1\mu m$.

The concentration of enzymes is chosen as $25 mM = 25 \cdot 10^{-3}M$

The specified values are to be given in SI units.


## Cases for which we know the analytical solution

**No reaction**: concentration equal to external concentration (trivial; convergence condition based on comparing reaction flux vs boundary flux not viable)

**Simple case of decay of compound X into compound Y, with no inner boundaries**: 

The reaction-diffusion equation in spherical coordinates is given by
$$ D \nabla^2 c - kc = D \frac{1}{r^2} \frac{\partial}{\partial r} \left( r^2 \frac{\partial c}{\partial r}\right) - kc = 0 $$
It is solved by
$$ c(r) = \frac{A\mathrm{exp}(-\lambda r) + B\mathrm{exp}(\lambda r)}{r}$$
with $\lambda = \sqrt{k/D}$.

The condition for internal reflexion $c^\prime(0) = 0$ leads to $A = -B$.

The condition at the outer membrane $c^\prime(R) = p/D \cdot (c_\mathrm{ext} - c(R))$ means 
$$ A = -\frac{\dfrac{p}{D}\,c_{\mathrm{ext}}\,R^2}
{e^{\lambda R}\left(\lambda R - 1 + \dfrac{pR}{D}\right)
+ e^{-\lambda R}\left(\lambda R + 1 - \dfrac{pR}{D}\right)} $$

**Simple case of decay of compound X into compound Y, with one inner boundary at $R^*$**:
Both individual segments (aka $[0, R^*)$ and $(R^*, R]$) follow
$$ c(r) = \frac{A\mathrm{exp}(-\lambda r) + B\mathrm{exp}(\lambda r)}{r}$$
with $\lambda = \sqrt{k/D}$.

We know
$$ c^\prime(0) = 0$$
$$ c^\prime(R^*_-) = c^\prime(R^*_+) = p / D \cdot (c(R^*_+)- c(R^*_-))$$
$$ c^\prime(R) = p/D \cdot (c_\mathrm{ext} - c(R))$$

We define 
$$ c(r)=
\begin{cases}
c_1(r), & 0\le r<R^*\\[6pt]
c_2(r), & R^*<r\le R
\end{cases}$$

Similarly to above, reflexion at $r=0$ leads to
$$ c_1(r)=\frac{S\,\sinh(\lambda r)}{r} $$

$$ c_2(r)=\frac{A e^{-\lambda r}+B e^{\lambda r}}{r} $$

We define
$$ s=\sinh(\lambda R^*), 
\qquad
c=\cosh(\lambda R^*), 
\qquad
\beta=\frac{p}{D},
\qquad
\alpha=\frac{pR}{D}.
$$
The final expression is given by
$$ B=\rho A $$

$$ \rho
=
e^{-2\lambda R^*}
\frac{
D\left(R^{*2}\lambda^2 c+R^*\lambda(c-s)-s\right)
+pR^{*2}\lambda(s+c)
}{
D\left(R^{*2}\lambda^2 c-R^*\lambda(s+c)+s\right)
+pR^{*2}\lambda(s-c)
}. $$

$$ A
=
\frac{\beta c_{\mathrm{ext}} R^2}{
e^{-\lambda R}(\alpha-\lambda R-1)
+
\rho e^{\lambda R}(\alpha+\lambda R-1)
},
$$


## How to run
Best to use conda:
1. make a new environment e.g. named snakemake_env
    ```bash
    conda create -n snakemake_env -c conda-forge -c bioconda snakemake
    ```
2. activate the new environment
    ```bash
    conda activate snakemake_env
    ```
3. follow the instructions in the Snakefile to run the scripts


## TODO:
- add ratios end_to_start_species
- add comparison to analytic solutions

## README Under construction

$$
 $$
is below some predefined threshold (e.g. 0.01).
reaction_flux = calculate_reaction_term(current_species_concentrations, region, n, species)
        reaction_fluxes[species] += 4 * np.pi * reaction_flux * r**2 * DELTA_R
    # Second, calculate flux from boundary with exterior
    # the flux is positive if the concentration on the exterior is larger than on the interior at r=R
    boundary_fluxes = {species: 0
        for species in REACTION_NETWORK.species}
    for species in REACTION_NETWORK.species:
        #print(species.name, species.permeability_constant, species.external_concentration, current_species_concentrations[NUM_REGIONS-1][NUM_MESH_POINTS_IN_REGIONS[NUM_REGIONS-1]-1][species])
        boundary_fluxes[species] = 4 * np.pi * R**2 * species.permeability_constant * (
            species.external_concentration
            - current_species_concentrations[NUM_REGIONS-1][NUM_MESH_POINTS_IN_REGIONS[NUM_REGIONS-1]-1][species])

relative_deviation = abs(reaction_fluxes[species] + boundary_fluxes[species]) / max(abs(boundary_fluxes[species]), abs(reaction_fluxes[species]))
