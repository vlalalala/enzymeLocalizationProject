# Numerical solver of chemical reactions within concentric spherical semi-permeable membranes

This code is for calculating the steady state concentrations of substances placed in a spherically symmetrical system of semi-permeable membranes, where the concentration of  each of substances outside the exterior membrane is kept constant. Spontaneous and enzymatic reactions can be defined (enzymes can be placed in the regions between semi-permeable membranes) and an arbitrary number of semi-permeable membranes can be used. The boundary problem is solved numerically through the Newton method. The solver comes to an end once it is not possible to further reduce the norm of the residual with the given mesh or once the net flux (reaction flux + boundary flux) of each species is "close enough" to 0 (see below for more details).

<img src="examples/simple_decay_with_two_inner_boundaries/combined_000001/newton_iterations.gif" loop=infinite>

## Summary

The governing equations for each species $q$ are reaction-diffusion equations of the form
```math
\frac{\partial q}{\partial t} = D_q \nabla^2 q + R_q
```
with $D_q$ the diffusion constant and $R_q$ the reaction term resulting from chemical reactions between the different species and/or enzymes,
and with the flux at the boundaries given by
```math
J = p \cdot (q(r^-) - q(r^+))
```
where $q(r^-)$ and $q(r^+)$ are the concentrations of $q$ at either side of the boundary. For each species, $q(R^+)$ is constant.

The steady state distribution of $q$ is computed numerically.

The interval $[0, R]$ is discretized by setting equally spaced mesh points. The position of each inner membrane is set to the closest mesh point position. For $N$ inner membranes defined, there are $N+1$ regions. Each region is defined by the mesh points within its bounds, including those at the bounding membranes. (Therefore, at the mesh positions where a inner boundary is at, there are in actuality 2 mesh points.)

The solution is assumed to have converged well enough if the net flux computed (the reaction flux minus the flux through the boundary) is close to zero for each species (i.e. mass conservation for each species, see theory part).

This work continues previous work by Hinzpeter et al: Optimal Compartmentalization Strategies for Metabolic Microcompartments, by Hinzpeter et al. (Biophysical Journal, 2017)

## How to use the code
The code uses `snakemake` for system management in order to be able to run reproducible and scalable data analyses. 
It is easy to use in combination with Anaconda/Miniconda. For usage with Miniconda:

1. make a new environment (e.g. named `snakemake_env`) and install snakemake
    ```bash
    conda create -n snakemake_env -c conda-forge -c bioconda snakemake
    ```
2. activate the new environment
    ```bash
    conda activate snakemake_env
    ```
3. follow the instructions in the Snakefile to run the scripts

## Parameters in the model
The parameters in the model are those defined in the following files (with example values introduced):

```yaml
# species.csv
```
| name    | diffusion_constant | external_concentration | permeability_constant |
| -------- | ----------------- | ---------------------- | --------------------- | 
| X  |  6.6e-10   | 90e-8 | 25e-6|
| Y  |  8.2e-10   | 0 | 15e-6|
| Z  |  4.2e-10   | 0 | 5e-6| 


```yaml
# spontaneous_reactions.csv
```
| start_species    | end_species | ratio_endtostart_species | k |
| -------- | ----------------- | ---------------------- | --------------------- | 
|  X |  Y   | 1:1| 1e-3| 


```yaml
# enzymes.csv
```
| name    | quantity | regions |
| -------- | ------- | -------- |
|  A  | 25e-2| [0,1] |

(A given quantity of) enzymes are placed within the regions specified. The concentration (in volume) is uniform.
```yaml
# enzymatic_reactions.csv
```
| start_species    | end_species | ratio_endtostart_species | enzyme | k_cat | k_M | hill |
| -------- | ------- | -------- | -------- | ------- | -------- | -------- |
| X |  Z   | 1:1| A  |  0.75   | 1.25e-4 | 1 |

Enzymatic reaction are defined given Michaelis-Menten kinetics. Hill exponents are added for cooperativity.

Automatic checks are included such that it is ensured that all species and enzymes that take place in reactions have their properties defined.

```yaml
# parameters_geometry.yaml
geometry_config:
    internal_membrane_relative_radii: [ 0.3, 0.7 ]
    outer_membrane_radius: 1.0e-5
```

```yaml
# parameters_solver_input.yaml

geometry_parameters: 
    num_mesh_points: 25

newton_parameters:
    override_adaptive_method: false

adaptive_step_parameters: 
    initial_alpha: 1.0
    alpha_min: 1.0e-3
    alpha_max: 10
    gamma_inc: 1.15
    gamma_dec: 0.5
    max_num_accepted_successive_unsuccessful_steps: 10

output_options:
    log_convergence_progress: true
    save_data_every: 100
    create_gif_with_saved_data: true
    log_iteration_info_every: 100
    delete_data_at_the_end: true
    plot_iteration_data_during_simulation: false

variables_to_save:
    save_F_vector: false
    save_F_vector_norm: false
    save_J_matrix: false
    save_du_vector: false
    save_concentrations: true
    save_du_vector_max: false
```

```yaml
# parameters_solver_params.yaml 

convergence_parameters:
    tol_relative_value: 1
    tol_absolute_factor: 1
    tol_residual_factor: 1
    tol_relative_flux_deviation: 0.01

newton_parameters:
    check_convergence_every: 100
```

The source code defined in the snakemake rules runs on files with this format.

**Defining the phase space spanned by combinations of parameter values:**

In order to efficiently run simulations testing out different regions in phase space, it is possible to provide a set of values for each of the parameters in the model. Individual simulations can then be run with each combination from the cartesian product of all these sets. 

To do so:
1. Create copies of the template files in the `/src` folder by running `python src/_create_parameters_template.py path_to_new_folder`.
2. Modify the entries in the `.yaml` and `.csv` files such that each entry contains a list of all the values that entry must take. (If all simulations share one same value for a given parameter, the list has length 1).
3. Create the files with the same format as above by running `python src/_create_templates_expanded.py path_to_new_folder` and `python src/_create_phase_space.py path_to_new_folder`. The different combinations are in folders that start with `combination*`.

## Theory
In spherical symmetry (assuming no angular dependence), the steady-state reaction-diffusion equation becomes

```math
\frac{1}{r^2} \frac{d}{dr} \left( r^2 \frac{du}{dr} \right) + \frac{1}{D_u} R_u(u,v...) = 0 \ ,
```
which simplifies to
```math
\frac{d^2u}{dr^2} + \frac{2}{r} \frac{du}{dr} + \frac{1}{D_u} R_u(u,v...) = 0 \ .
```

(Note the singularity at $r=0$).

Using
```math
\begin{cases}
y_0 = u \\
y_1 = \frac{du}{dr} \\\end{cases}
```
we define each of the $2^\mathrm{nd}$ order differential equations (one for each species) as 2 $1^\mathrm{st}$ order differential equations of the form
```math
\begin{cases}
\frac{dy_0}{dr} = y_1 \\
\frac{dy_1}{dr} = -\frac{2}{r} y_1 - \frac{1}{D} R_y(y_0, z_0...)
\end{cases}
```

As for the boundary conditions:

At $r=0$ we use reflection
```math
\frac{du}{dr}(0) = 0
```
such that
```math
\left(\frac{dy_0}{dr}(0) = \right)\ \  y_1(0) = 0
```
and at $r=R$ we have
```math
\frac{du}{dr}(R) = [u_\mathrm{ext} − u(R)] \cdot \frac{p_u}{D_u}
```
such that
```math
y_1(R) = (u_{ext} - y_0(R)) \cdot \frac{p_u}{D_u}
```
with $p_u$ and $D_u$ the permeability to each membrane and diffusion constant for substance $u$, respectively.

**Generally**:

Given $m$ species and $n$ compartments (i.e. $n-1$ inner membranes at $r_i$ with $i \in \{1,2,...,n-1\} $) we have $m \cdot [2 + (n-1)\cdot 2] = 2mn$ boundary conditions (since for each species we have one BC at $r=0$, one BC at $r=R$, one BC at each side of each inner membrane).

In order to use the Newton method to solve the whole system, we define $2\cdot m\cdot n$ variables, with $k\in \{0, 1, ..., n-1\}$

```math
v_k^{(m)} = u_k^{(m)}
```
```math
w_k^{(m)} = \frac{\mathrm{d}u_k^{(m)}}{\mathrm{d}r}
```

with 

```math
\frac{\mathrm{d}v_k^{(m)}}{\mathrm{d}r} = w_k^{(m)}
```
```math
\frac{\mathrm{d}w_k^{(m)}}{\mathrm{d}r} = -\frac{2}{r} w_k^{(m)}-\frac{1}{D}R(\vec{v}_k)
```

Boundary conditions are given by
```math
w^{(m)}(0) = 0
```
(Neumann no-flux BC)
```math
w^{(m)}(R) = \left(u_\mathrm{ext} - v^{(m)}(R)\right) \cdot \frac{p^{(m)}}{D}
```
(Robin BC)
```math
w_i^{(m)}(r_i^+) =  \frac{p^{(m)}}{D} \cdot \left( v_{i}^{(m)}(r_i^+) - v_{i-1}^{(m)}(r_i^-)\right)
```
```math
w_{i-1}^{(m)}(r_i^-) =  \frac{p^{(m)}}{D} \cdot \left( v_{i}^{(m)}(r_i^+) - v_{i-1}^{(m)}(r_i^-)\right)
```

The concentration is not continuous at the membranes, but the flux is continuous!

We used

```math
f'(x_i) = \frac{1}{2h}(f_{i+1}-f_{i-1}) + O(h^2)
```
```math
f''(x_i) \approx \frac{1}{h^2}(f_{i+1}-2f_i + 
f_{i-1})
```
for the discretization.

**Given multiple regions**:

In each region $k$ we have 

```math
\frac{1}{r^2} \frac{d}{dr} \left( r^2 \frac{du_k}{dr} \right) + \frac{1}{D_k} R(u) = 0 ,
```

```math
D_k \, \frac{\mathrm{d}u_k}{\mathrm{d}r} \big|_{r = R_k^-}
= D_{k+1} \, \frac{\mathrm{d}u_{k+1}}{\mathrm{d}r} \big|_{r = R_k^+}
= p_k \left( u_{k+1}(R_k^+) - u_k(R_k^-) \right)
```


**Fluxes**

**Finding 
To find the total number of moles $n$ within the volume with radius R, calculate
```math
n = \int_0^R C(r) \cdot 4\pi r^2 \mathrm{d}r
```

Reminder that to convert $C(r)$ from $M$ (moles per litre) to moles per $m^{3}$: $M = \frac{moles}{m^3}\cdot 10^{-3}$.

So 
```math
n = \int_0^R C(r) \cdot 10^{3}\cdot 4\pi r^2 \mathrm{d}r
```
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
```math
D \nabla^2 c - kc = D \frac{1}{r^2} \frac{\partial}{\partial r} \left( r^2 \frac{\partial c}{\partial r}\right) - kc = 0
```

It is solved by
```math
c(r) = \frac{A\mathrm{exp}(-\lambda r) + B\mathrm{exp}(\lambda r)}{r}
```
with $\lambda = \sqrt{k/D}$.

The condition for internal reflexion $c^\prime(0) = 0$ leads to $A = -B$.

The condition at the outer membrane $c^\prime(R) = p/D \cdot (c_\mathrm{ext} - c(R))$ means 
```math
A = -\frac{\dfrac{p}{D}\,c_{\mathrm{ext}}\,R^2}
{e^{\lambda R}\left(\lambda R - 1 + \dfrac{pR}{D}\right)
+ e^{-\lambda R}\left(\lambda R + 1 - \dfrac{pR}{D}\right)}
```

**Simple case of decay of compound X into compound Y, with one inner boundary at $R^*$**:
Both individual segments (aka $[0, R^*)$ and $(R^*, R]$) follow
```math
c(r) = \frac{A\mathrm{exp}(-\lambda r) + B\mathrm{exp}(\lambda r)}{r}
```
with $\lambda = \sqrt{k/D}$.

We know
```math
c^\prime(0) = 0
```
```math
c^\prime(R^*_-) = c^\prime(R^*_+) = p / D \cdot (c(R^*_+)- c(R^*_-))
```
```math
c^\prime(R) = p/D \cdot (c_\mathrm{ext} - c(R))
```

We define 
```math
c(r)=
\begin{cases}
c_1(r), & 0\le r<R^*\\[6pt]
c_2(r), & R^*<r\le R
\end{cases}
```


Similarly to above, reflexion at $r=0$ leads to
```math
c_1(r)=\frac{S\,\sinh(\lambda r)}{r}
```
```math
c_2(r)=\frac{A e^{-\lambda r}+B e^{\lambda r}}{r}
```

We define
```math
s=\sinh(\lambda R^*), 
\qquad
c=\cosh(\lambda R^*), 
\qquad
\beta=\frac{p}{D},
\qquad
\alpha=\frac{pR}{D}.
```
The final expression is given by
```math
B=\rho A 
```

```math
\rho
=
e^{-2\lambda R^*}
\frac{
D\left(R^{*2}\lambda^2 c+R^*\lambda(c-s)-s\right)
+pR^{*2}\lambda(s+c)
}{
D\left(R^{*2}\lambda^2 c-R^*\lambda(s+c)+s\right)
+pR^{*2}\lambda(s-c)
}.
```


```math
A
=
\frac{\beta c_{\mathrm{ext}} R^2}{
e^{-\lambda R}(\alpha-\lambda R-1)
+
\rho e^{\lambda R}(\alpha+\lambda R-1)
},
```





## TODO:
- add ratios end_to_start_species
- add comparison to analytic solutions

## README Under construction


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
