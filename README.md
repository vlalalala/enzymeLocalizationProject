---
fontsize: 10pt
papersize: a4
geometry: margin=0.75in
---

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

**IMPORTANT**: ratios different to 1:1 between reactants and products is not yet implemented!

## Optimization
We optimize the flux of some substance by allocating the different enzymes across the different regions.

The number of free parameters is the following: (# enzyme types -1) + (# enzymes) x (# regions - 1) + (# regions - 1): 
    1. The first factor represents the degrees of freedom we get from allocating a total fixed quantity of enzyme towards the different enzymes.
    2. The second factor is the degrees of freedom from each enzyme allocation to the different regions has to sum up to 100%.
    2. The third (# regions - 1) comes from having to find the optimal position of the (# regions - 1) number of inner membranes.

Softmax function:
```math
x_i = \frac{e^{z_i}}{\sum_j e^{z_j}}
```

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

**Running on a SLURM cluster** (not tested yet):
1. make an environment/install the following packages
    ```bash
    conda install -c conda-forge -c bioconda snakemake snakemake-executor-plugin-slurm
    ```
2. Run once from the command line:
    ```bash
    snakemake \
    --executor slurm \
    --jobs 100 \
    --use-conda \
    --default-resources \
    --rerun-incomplete
    ```
    or using profile data, run using `snakemake --profile profiles/slurm`
Alternatively to step 2, if ssh connection can break, submit through `sbatch run_all.sh`



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
    internal_membrane_relative_radii: [ 0.3, 0.7 ] # the positions of inner membranes in terms of the outer-most radius. list
    outer_membrane_radius: 1.0e-5 # given in meters
```

```yaml
# parameters_solver_input.yaml

geometry_parameters: 
    num_mesh_points: 25 # number of different mesh positions
                        # (the total number of mesh points is larger due to
                        # duplicate mesh points at inner boundaries)

newton_parameters:
    override_adaptive_method: false # an adaptive method for the step in concentrations
                                    # is used as the default

adaptive_step_parameters: 
    initial_alpha: 1.0 # alpha is the scaling factor for the computed step 
                        # in concentrations via the Newton method
    alpha_min: 1.0e-3
    alpha_max: 10
    gamma_inc: 1.15
    gamma_dec: 0.5 # factor by which alpha is multiplied when the residual norm does not decrease
    max_num_accepted_successive_unsuccessful_steps: 10 # used to stop iterations of Newton once the
        # most accurate solution with the given number of mesh points is found

output_options:
    log_convergence_progress: true # if true, logs relative net flux for each species
    save_data_every: 100 # save quantities specified by variables_to_save given 
                            # this number of iterations; if nothing should be saved, set to 0
    create_gif_with_saved_data: true # to create a gif with the saved concentrations
                                     # at the end of the simulation
    log_iteration_info_every: 100 # log inforation about (adaptive) step size and 
                                    # number of successive unsuccessful steps
    delete_data_at_the_end: true # if true, deletes files specified by variables_to_save at the end
    plot_iteration_data_during_simulation: false # if true, plots the .png file for the latest
                                                # saved iteration with the latest calculated concentration

variables_to_save:
    save_F_vector: false # save residual vector as a .txt
    save_F_vector_norm: false # save norm of residual vector as a .txt
    save_J_matrix: false # save the (sparse) jacobian matrix as a .txt
                            # (saves only the non-zero elements with their value and position)
    save_du_vector: false # save the vector that encodes the change in the concentrations
    save_concentrations: true # save the concentrations (needed for creating the gif)
    save_du_vector_max: false # save the maximum step change in the concentrations 
```

```yaml
# parameters_solver_params.yaml 

convergence_parameters:
    tol_relative_value: 1 # factor by which the initial norm of the residual must be reduced 
    tol_absolute_factor: 1 # obsolete, similar to above; read code
    tol_residual_factor: 1 # obsolete, similar to above; read code
    tol_relative_flux_deviation: 0.01 # upper threshold for the relative net flux
                                        # with which the balance condition is accepted

newton_parameters:
    check_convergence_every: 100 # every how many iterations to check the convergence coditions
```

The source code defined in the snakemake rules runs on files with this format.

**Defining the phase space spanned by combinations of parameter values:**

In order to efficiently run simulations testing out different regions in phase space, it is possible to provide a set of values for each of the parameters in the model. Individual simulations can then be run with each combination from the cartesian product of all these sets. 

To do so:
1. Create copies of the template files in the `/src` folder by running `python src/_create_parameters_template.py path_to_new_folder`.
2. Modify the entries in the `.yaml` and `.csv` files such that each entry contains a list of all the values that entry must take. (If all simulations share one same value for a given parameter, the list has length 1).
3. Create the files with the same format as above by running `python src/_create_templates_expanded.py path_to_new_folder` and `python src/_create_phase_space.py path_to_new_folder`. The different combinations are in folders that start with `combination*`.

IMPORTANT: it is very important that the ratios between end and start species are written in `" " ` in order that it is read as a string and not a sexagesimal number.

## Theory
We start from the reaction-diffusion system in a steady state
```math
0 = D_s \nabla^2 c_s + R_s(\bm{c})
```
In spherical symmetry (assuming no angular dependence), the steady-state reaction-diffusion equation (for each species) becomes

```math
\frac{1}{r^2} \frac{\mathrm{d}}{\mathrm{d}r} \left( r^2 \frac{\mathrm{d}c_s}{\mathrm{d}r} \right) + \frac{1}{D_s} R_s(\bm{c}) = 0 \ ,
```
which simplifies to
```math
\frac{\mathrm{d}^2 c_s}{\mathrm{d}r^2} + \frac{2}{r} \frac{\mathrm{d}c_s}{\mathrm{d}r} + \frac{1}{D_u} R_s(\bm{c}) = 0 
```
(note the singularity at $r=0$).

Since there are concentration jumps at the membranes, we discretize the interval $[0,R]$ with the following procedure:
<div style="background-color: white; display: inline-block; padding: 10px;">
  <img src="notes/discretization.png">
</div>

We have $M = R/h+1$ "positions" for mesh points.
Given $S$ species and $N$ compartments (i.e. $N-1$ inner membranes) we have a system of $(M+N+1) \cdot S$ equations.

There are $S \cdot [2 + (N-1)\cdot 2] = 2SN$ boundary conditions (since for each of the S species we have one BC at $r=0$, one BC at $r=R$ and one BC at each side of each inner membrane).

For the discretization of space, we use
```math
f^\prime(x_i) = \frac{1}{2h}(f_{i+1}-f_{i-1}) + O(h^2) \\[1em]
f^{\prime\prime}(x_i) \approx \frac{1}{h^2}(f_{i+1}-2f_i + f_{i-1})
```
with $h$ the distance between adjacent mesh points (within the same region).

**Interior points**

For each interior point, we have
```math
0 = D_s \cdot \left( \frac{c_s^\mathrm{right} - 2 \cdot c_s^\mathrm{center} + c_s^\mathrm{left}}{h^2} + \frac{c_s^\mathrm{right} - c_s^\mathrm{left}}{h \cdot r}\right) + R_s\left(\bm{c}^\mathrm{center}\right)
```
where $r$ is the distance from the origin to the raidus of the interior point.

**Treatment of points at $r=0$**

For the points at $r=0$, we expand $c(r)$ in a Taylor series around $0$:
```math
c(r) = c(0) + c^\prime (0) r + \frac{1}{2} c^{\prime\prime}r^2 + \mathcal{O}(r^3)
```
s.t.
```math
c^\prime(r) = c^{\prime\prime}(0)r + \mathcal{O}(r^2)
```
and
```math
c^{\prime\prime}(r) = c^{\prime\prime}(0) + \mathcal{O}(r)
```

We have
```math
c^\prime (0) = 0
```
due to symmetry.

Recall
```math
\frac{\mathrm{d}^2 c_s}{\mathrm{d}r^2} + \frac{2}{r} \frac{\mathrm{d}c_s}{\mathrm{d}r} + \frac{1}{D_u} R_s(\bm{c}) = 0 
```
The second term can be written in the following way:
```math
\frac{2}{r}c^\prime(r) = \frac{2}{r}(c^{\prime\prime}(0)r + \mathcal{O}(r^2)) = 2c^{\prime\prime}(0)+\mathcal{O}(r)
```
s.t.
```math
\left(\frac{\mathrm{d}^2 c_s}{\mathrm{d}r^2} + \frac{2}{r} \frac{\mathrm{d}c_s}{\mathrm{d}r}\right) \Bigr|_{r=0} \approx 3 c^{\prime\prime}(0)
```




```math
\nabla^2 c(0) = 3 \frac{\partial^2 c}{\partial r^2}(0)
```
(It is also possible to derive this using l'Hospital rule.)
Using $c_{-1} = c_1$ we get
```math
\nabla^2 c_s(0) = 3 \cdot D_s / h^2 \cdot 2 \cdot (c_s^\mathrm{right} - c_s^\mathrm{center})
```
such that
```math
0 = 3 \cdot D_s / h^2 \cdot 2 \cdot (c_s^\mathrm{right} - c_s^\mathrm{center}) + R_s\left(\bm{c}^\mathrm{center}\right)
```

**Points at the semipermeable membranes**

For the points at the right of the semipermeable membrane (i.e. the left-most points of the region to the right of the membrane) we use
```math
0 = D_s \cdot \frac{c_s^\mathrm{right} - c_s^\mathrm{center,+}}{h} - p_s \cdot (c_s^\mathrm{center, +} - c_s^\mathrm{center, -})
```
where $c_s^\mathrm{center,+}$ is the concentration at the membrane on its right side and $c_s^\mathrm{center,-}$ is the concentration on the left side.

For the points at the left of the semipermeable membrane (i.e. the right-most points of the region to the left of the membrane), we use
```math
0 = D_s \cdot \frac{c_s^\mathrm{center, -} - c_s^\mathrm{left}}{h} - p_s \cdot (c_s^\mathrm{center, +} - c_s^\mathrm{center, -})
```
analogously.

Question: should I also put reaction flux here?

This gives a nonlinear algebraic system
$F(\vec(x)) = 0$ where $\vec(x)$ contains all species concentrations at all grid points.

The Jacobian $J$ is given by
```math
J_{ij} = \frac{\partial F_i}{\partial x_j}
```


**Using the Newton method**

Normally
```math
J(\bm{x}^{(k)}) \delta \bm{x}^{(k)} = -F(\bm{x}^{(k)}) \\

\bm{x}^{(k+1)} = \bm{x}^{(k)} + \delta \bm{x}^{(k)}
```
or (equivalently)
```math
\delta \bm{x} = - J^{-1} F
```
(since J is invertible, multiplying both sides of the first equation by $J^(-1)$ through the left)

However: Newton method only works (in this way) close to the solution.

Therefore, we use an adaptive step size:
Following the paper *An adaptive Newton-method based on a dynamical systems approach* by Amrein and Wihler, we use

```math
x_{n+1} = x_n - t_n J(x_n)^{-1}F(x_n)\\[1em]
t_n = \mathrm{min}\left(\sqrt{\frac{2\tau}{||\mathrm{N}_\mathrm{F}(x_n)||_X}} , 1\right)
```
where $\mathrm{N}_\mathrm{F}(x_n) = -J^{-1}F(x)$ and $\tau$ is some tolerance. We define $\tau$ dynamically


```math
\begin{align*}
&\text{Start: } u, \tau, \|F\|_{\text{last}} \\
&\text{Compute residual and Jacobian: } F(u), J(u) \\
&\Delta u_N = - J(u)^{-1} F(u), \quad \|\Delta u_N\| = \text{norm}(\Delta u_N) \\
&t_n = \min\Big(\sqrt{\frac{2 \tau}{\|\Delta u_N\|}}, 1\Big), \quad \Delta u = t_n \Delta u_N \\
&u_{\text{trial}} = u + \Delta u \\
&\text{Check positivity:} \\
&\quad \text{If any } u_{\text{trial}} < 0: \\
&\qquad \tau_{\text{new}} = \tau \cdot \gamma_{\text{dec}} \\
&\qquad \text{If } \tau_{\text{new}} < \tau_{\min} \text{ → ERROR (negative values)} \\
&\qquad \text{Return } (u, \tau_{\text{new}}, \|F\|_{\text{last}}) \\
&\text{Compute new residual: } F(u_{\text{trial}}), \quad \|F\|_{\text{new}} = \text{norm}(F(u_{\text{trial}})) \\
&\text{Decision:} \\
&\quad \text{If } \|F\|_{\text{new}} < \|F\|_{\text{last}} \quad \text{(successful step)} \\
&\qquad \tau_{\text{new}} = \min(\tau_{\max}, \tau \cdot \gamma_{\text{inc}}) \\
&\qquad \text{Return } (u_{\text{trial}}, \tau_{\text{new}}, \|F\|_{\text{new}}) \\
&\quad \text{Else (unsuccessful step)} \\
&\qquad \tau_{\text{new}} = \tau \cdot \gamma_{\text{dec}} \\
&\qquad \text{If } \tau_{\text{new}} < \tau_{\min} \text{ → ERROR (cannot reduce residual further)} \\
&\qquad \text{Return } (u, \tau_{\text{new}}, \|F\|_{\text{new}})
\end{align*}

```

**Calculating the net flux**

As a convergence condition, we establish that the (total) net flux should be close to 0.

We define the flux coming through the exterior membrane to be
```math
\Phi_\mathrm{ext}^{(m)} = 4 \cdot \pi \cdot p_u\cdot (v_\mathrm{ext}^{(m)} - v^{(m)}(R))
```
such that it is positive if there is a flux into the sphere and negative if there is a flux towards the outside of the sphere.

We estimate the reaction flux for each region by calculating the reaction fluxes caused by the local concentrations at each mesh point (including the mesh points at either end of the interval), interpolating the reaction fluxes at the positions between each of these mesh points and adding these reaction fluxes, weighted by the volume of the hollow sphere spanned by the two elements from which the interpolation was calculated.
```math
\Phi_\mathrm{react}^{(m)} = \sum_{i=0}^{N-1} \frac{R^{(m)}_i + R^{(m)}_{i+1}} {2}\cdot \frac{4}{3} \cdot \pi \cdot \left(r_{i+1}^3 - r_i^3\right)
```
with the index $i$ traversing through the mesh points within a region with increasing $r$ and $R_i^{(m)}$ the computed reaction term from the local concentrations of all species and enzymes.
The total reaction flux is the sum from the reaction fluxes within each region.

We define a relative net flux $\Phi^{(m)}$ for each species $m$ through
```math
\Phi^{(m)} = \frac{|\Phi_\mathrm{ext}^{(m)} + \Phi_\mathrm{react}^{(m)}|}{\mathrm{max}(|\Phi_\mathrm{ext}^{(m)}|, |\Phi_\mathrm{react}^{(m)}|)}
```
Balance between the flux created through interaction with the outside and the flux created through reactions within the sphere must be ensured for the steady state. We consider that the steady state has been found numerically once this balance is smaller than some small $\epsilon > 0$.

It is important to note that the net flux may not cross the threshold given by $\epsilon$ if the step size between neighboring mesh points is too large.


## How to make an informed initial guess for the concentrations

### Plan A: Solve multiple simulations, each with higher reaction

For no reaction, the steady state concentration is equal to the external concentration for each species.

We use the Thiele Modulus $\phi = R \cdot \sqrt{k/D}$ to assess how fast reactions are with respect to diffusion ($\phi \ll 1$ means that diffusion is fast; $\phi \gg 1$ means that reactions are fast).

We run a number $N_\mathrm{sim}$ of simulations.
$k$ and $k_\mathrm{cat}$ modified.

For each simulation:
- For each species:
    - Clip $k$ and $k_\mathrm{cat}$ by the same factor s.t. the largest reaction timescale is ca. 10x of diffusion timescale -> will produce a result close to the external concentration, but slightly shifted. 

    The reaction (total) timescale relates to the timescales of the individual reactions is given by 
    ```math
    \frac{1}{\tau_\mathrm{reaction}} = \frac{1}{\tau_\mathrm{1}} + \frac{1}{\tau_\mathrm{2}}.
    ```

The following condition must hold:
```math
x \cdot \gamma_\mathrm_{inc} \cdot \gamma_\mathrm_{dec} > x
```
so that decreasing the factor does not bring us to a factor x that is below one that was successful.
We choose a stronger condition: If a simulation fails to converge, the new factor chosen must be
in the middle between the one that failed and the previous one ()
```math

```


### Plan B: Estimate initial value
#### Species gets consumed through spontaneous reaction or through enzyme in linear regime

Assume a simple reaction of 1st order
```math
D_s\cdot\frac{1}{r^2} \frac{\mathrm{d}}{\mathrm{d}r} \left( r^2 \frac{\mathrm{d}c_s}{\mathrm{d}r} \right) - k \cdot c_s = 0 \ ,
```
Using
```math
u(r) = r \cdot c_s(r)
```
we get
```math
D_s \cdot \frac{\mathrm{d}^2 u}{\mathrm{d}r^2}- k \cdot u = 0
```
s.t. the general solution is
```math
u(r) = A \cdot \mathrm{sinh}(\phi \cdot r/R) + B \cdot \mathrm{cosh}(\phi \cdot r/R)
```
where $\phi = R \cdot \sqrt{k/D}$ is the Thiele modulus.

At $r=0$ the concentration must be finite. This means that $u(0)$ must be equal to zero (s.t. $c(0) = u(0)/0$ does not blow up). Because $\mathrm{cosh}(0) = 1$, we need $B=0$.


We use $c(R) = c_\mathrm{ext}$ (which doesn't have to be the case). This means that $u(R) = A \cdot \mathrm{sinh}(\phi) = c_\mathrm{ext} / R$. Thus, $A = c_\mathrm{ext} \cdot R / \mathrm{sinh}(\phi)$.

The solution thus takes the form
```math
c(r) = c_\mathrm{ext} \cdot (R/r) \cdot \frac{\mathrm{sinh}(\phi \cdot r/R)}{\mathrm{sinh}(\phi)}
```

$k$ is the net rate with which the substance is consumed, either through a spontaneous reaction of 1st order or in the linear regime of an enzymatic reaction ($S\ll k_M$; using $k_\mathrm{eff} = k_\mathrm{cat}\cdot E / k_M$).

What happens if enzyme is saturated ($S\gg k_M$)?

#### Species gets consumed through spontaneous reaction, through enzyme in linear regime and through enzyme in saturated regime
For this, we solve

```math
D_s\cdot\frac{1}{r^2} \frac{\mathrm{d}}{\mathrm{d}r} \left( r^2 \frac{\mathrm{d}c_s}{\mathrm{d}r} \right) - k \cdot c_s -k_\mathrm{cat}\cdot E = 0 \ ,
```
with $E$ the enzyme concentration. The solution then is
```math
c(r) = (c_\mathrm{ext} + k_\mathrm{cat} \cdot E / k) \cdot (R/r) \cdot \mathrm{sinh}(\phi \cdot r/R) / \mathrm{sinh}(\phi)  -  k_\mathrm{cat} \cdot E / k
```
where we are careful to not count $k_\mathrm{cat}$ more than once (i.e. do not consider it into $\phi$).

**Important**: I    t might make sense to use the linear regime approximation regardless with the external concentration. If we are (in reality) nevertheless in the other regime, it will overestimate the depletion, though.


#### Species is produced
For the initial guess, the product of the mirror of the substrate.

## Literature on numerics used:
- C.T.Kelley: Iterative Methods for Linear and Nonlinear Equations $\rightarrow$ "condition number" 
- W.H.Press, S.A.Teukolsky, W.T. Vetterling, B.P.Flannery: Numerical Recipes in C, The Art of Scientific Computing $\rightarrow$ book recommended by Uli
- L.N.Trefethen, D. Bau: Numerical Linear Algebra (lectures) $\rightarrow$ "condition number" (Lecture 12), has the analysis on perturbation.  


## Convergence
### Criterion 1
A common convergence criterion (see Lemma 1.1.1 from C.T.Kelley) is the following:
When solving the problem $Ax=b$, we can use
```math
\frac{\lVert r_k\rVert}{\lVert b \rVert} < \tau
```
where $r_k$ is the residual for the $k$-th iteration.

This begs the question on how big $\tau$ should be.

### Criterion 2
One can decide to stop iterating once the change in concentrations (for each concentration, element-wise) is below some threshold (tolerance):
```math
\max_{i} \left|\frac{\delta u_i}{u_i}\right| < \mathrm{tol}
```
Machine precision $\epsilon$ (also called machine epsilon) leads us to the following two points to consider:

1. For 64-bit floats, $\epsilon \approx 2.2\cdot 10^{-16}$. If $\delta u_i < u_i \cdot \epsilon$ for some $i$, $u_i$ will effectively not be updated. 
2. Additionally, one has to consider how accurately one can solve systems of equations given such a machine precision. 
(see Theorems 12.1 and 12.2 Trefethen): "If a problem $Ax=b$ contains an ill-conditioned matrix $A$, one must always expect to "lose $\mathrm{log}_{10} \kappa(A)$ digits" in computing the solution", with $\kappa(A)$ being the condition number $\lVert A \rVert \lVert A^{-1}\rVert$. The floor is then $\kappa(J) \cdot \epsilon$ and the tolerance must be slightly above (since the floor is not being crossed.)

If $\kappa(A)$ is small, $A$ is said to be well-conditioned; if $\kappa(A)$ is large, $A$ is ill-conditioned. 

```math
\frac{\lVert \delta x \rVert}{\lVert x \rVert} \le \kappa(A) \cdot \frac{\lVert \delta b \rVert}{\lVert b \rVert}
```

See below for a more detailed proof:
```math
J(\textbf{u}) \cdot \delta \textbf{u} = F(\textbf{u})
\rightarrow
\delta \textbf{u} = J(\textbf{u})^{-1} F(\textbf{u})\\

J(\textbf{u}) \cdot \delta \textbf{u}^\prime = F(\textbf{u}) + \delta F(\textbf{u})
\rightarrow
\delta \textbf{u}^\prime = J(\textbf{u})^{-1} (F + \delta F) \\

\delta (\delta u) = \delta u^\prime - \delta u = J^{-1} \delta F \\

\lVert \delta(\delta u) \rVert \le \lVert J^{-1}\rVert \cdot \lVert \delta F \rVert \\

\lVert F \rVert \le \lVert J \rVert \cdot \lVert \delta u \rVert \\

Dividing: \\
\lVert \delta(\delta u)\rVert / \lVert \delta u \rVert \le \lVert J^{-1} \rVert \cdot \lVert J \rVert \cdot \lVert \delta F \rVert / \lVert F \rVert \\

```

The condition number is large if the Jacobian has very different scales accros its rows or columns.

so small that it is insignificant given the order of magnitude of the concentration. (This means that `u_n + alpha * du == u_n` element_wise for the step $n$ (this should apply for each element). e.g. written through `np.max(np.abs(du / u)) < tol`)
    - Careful if `u_i ≈ 0`

Questions: 
- Newton decrement


## Typical values
see 
https://bionumbers.hms.harvard.edu/search.aspx?trm=diffusion+coefficient

**Diffusion coefficient**
- Diffusion coefficient of l-Trp estimated as $D = 6.6\cdot 10^{-6} cm^{2} s^{-1} = 6.6\cdot 10^{-10} m^{2} s^{-1}$ (estimated from https://pmc.ncbi.nlm.nih.gov/articles/PMC16526/).
- Diffusion coefficient in "Optimal Compartmentalization Strategies for Metabolic Microcompartments" (Hinzpeter et al) estimated as $D=1000 \mu m^2s^{-1} = 1.0 \cdot 10^{-9} m^2 s^{-1}$

**Decay constants**
Between 10^{-6} and 10^{3}.

The decay rates can be calculated from the half-lives of substances
```math
k = \frac{2}{t_{1/2}}
```
- For fast-signaling molecules:
    Half-life: $10^{-3}-10^2 s$, k between 1e-2 and 1e3 1/s
    Check references given by ChatGPT:
    Stryer et al., Biochemistry (signal transduction kinetics chapters)
    Purves et al., Neuroscience (synaptic clearance kinetics)
    Berridge (2006), Cell Calcium signaling dynamics
- Enzyme produced metabolites
    Half-life 1s-10^4 s; k from 1e-4 to 1 1/s
    Check references given by ChatGPT:
    Berg, Tymoczko, Stryer — Biochemistry
    Heinrich & Schuster (1996), The Regulation of Cellular Systems
    Noe et al. (2013), metabolic flux modeling in cells
- Protein degradation
    Half-life 10min - 48h; k from 1e-6 to 1e-3 1/s
    Check references given by ChagGPT:
    Schwanhäusser et al., Nature (2011) — global protein half-lives in mammalian cells
    Eden et al., Nature (2011)
    Cambridge Systems Biology textbook (Klipp et al.)


**Maximal enzyme concentration**
- Maximal enzyme concentration in "Optimal Compartmentalization Strategies for Metabolic Microcompartments" (Hinzpeter et al) estimated as $25 mM = 25 \cdot 10^{-3}M = 25 mol/m^3$.
- According to Zotter et al ("Quantifying enzyme activity in living cells"): "we were able to determine catalytic constants in cells with enzyme concentrations ranging from 0.01 to 10 μM and substrate concentration ranging from 1 to 100 μM. These are common concentrations for enzymes and substrates in the cell (34)". (34) is "BioNumbers—the database of key numbers in molecular and cell biology" by Milo et al.

**Enzyme kinetic parameters**
- Characteristics of VioA following https://pubs.acs.org/doi/pdf/10.1021/acsomega.3c08233?ref=article_openPDF: $k_\mathrm{M} = 125 \mu M = 125\cdot 10^{-3}mM = 125 \cdot 10^{-3} mol/m^3$ and $k_\mathrm{cat} = 0.75 s^{-1}$
- Full range of Typical K_m values is between 1e-7 and 1e-1
- Most enzymes have a K_m between 1e-6 and 1e-3
- Kcat is usually between 1 and 1e4

**External concentration**
- External concentration following https://www.cell.com/biophysj/fulltext/S0006-3495(16)34263-1?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0006349516342631%3Fshowall%3Dtrue : 
l-Trp we're assuming $s_0 = 25 \mu M = 25\cdot 10^{-6}M = 25\cdot 10^{-3}mol/m^3$.
- Important note to self:
    - it might make sense to take the external concentration to be within the range of $(0.1 - 10) \cdot k_M$
- According to Zotter et al ("Quantifying enzyme activity in living cells"): "we were able to determine catalytic constants in cells with enzyme concentrations ranging from 0.01 to 10 μM and substrate concentration ranging from 1 to 100 μM. These are common concentrations for enzymes and substrates in the cell (34)". (34) is "BioNumbers—the database of key numbers in molecular and cell biology" by Milo et al.

**Permeability constant**
- Maximal enzyme concentration in "Optimal Compartmentalization Strategies for Metabolic Microcompartments" (Hinzpeter et al) estimated in the Supplementary Material. They use $p = 90 \mu M s^{-1} = 90\cdot 10^{-6}M s^{-1} = 90\cdot 10^{-3} mol s^{-1}$ for some and a fifth of that for others.

**Vesicle size**
- Radius of $1\mu m$.
(Note:Half the volume is in inner 79% of radius)

The specified values are all to be given in SI units (including concentrations)!

**Enzyme quantity**
keep in mind that $q_\mathrm{max} = c_\mathrm{max} * 4/3 * \pi * R**3$, which for $R = 1^{-6}m$ and $c_\mathrm{max}=25 \cdot 10^{-3}mol/m^3$ leads to $q_\mathrm{max} = 1.05\cdot 10^{-19}mol$

## Metabolic control theory
- Cells tend to tune K_m relative to substrate concentration

- Near-phisiology matching: K_m approx substrate concentration -> maximizes sensitivity to concentration changes
- Low K_m (enzyme saturated): "buffering behavior" (constant flux)
- High K_m (enzyme acts as a sensor or only active at high substrate load)

Books:
- Heinrich & Rapoport (1974)
- Fell (1997), Understanding the Control of Metabolism

## Cases for which we know the analytical solution

**No reaction**: concentration equal to external concentration (trivial; convergence condition based on comparing reaction flux vs boundary flux not viable)

```math
D \nabla^2 c(r) = D \frac{1}{r^2} \frac{\partial}{\partial r} \left( r^2 \frac{\partial c(r)}{\partial r}\right) = 0
```
If $c(r)=C$, we get $0=0$.


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


**Comparison with case given non-spherical, purely 1D diffusion**
Dropping the spherical symmetry, the steady-state reaction-diffusion equation on a line is
```math
D \frac{d^2 c}{d x^2} - k\cdot c = 0
```
The solution to that is
```math
c(x) = A\cosh\!\left(\frac{x}{\ell}\right) + B\sinh\!\left(\frac{x}{\ell}\right)
```
with $l = \sqrt{D/k}$.

$c'(0) = 0$ leads to
```math
    c'(x) = \frac{A}{\ell}\sinh\!\left(\frac{x}{\ell}\right) + \frac{B}{\ell}\cosh\!\left(\frac{x}{\ell}\right),
```
such that $B=0$.

$c'(R) = \frac{p}{D}\bigl(c_{\mathrm{ext}} - c(R)\bigr)$ leads to 

Substituting $c(x) = A\cosh(x/\ell)$:
```math
    \frac{A}{\ell}\sinh\!\left(\frac{R}{\ell}\right)
    = \frac{p}{D}\!\left[c_{\mathrm{ext}} - A\cosh\!\left(\frac{R}{\ell}\right)\right]
```
```math
A\left[\frac{1}{\ell}\sinh\!\left(\frac{R}{\ell}\right)
          + \frac{p}{D}\cosh\!\left(\frac{R}{\ell}\right)\right]
    = \frac{p}{D}\,c_{\mathrm{ext}},
```
which yields
```math
    A = \frac{\dfrac{p}{D}\,c_{\mathrm{ext}}}
             {\dfrac{1}{\ell}\sinh\!\left(\dfrac{R}{\ell}\right)
              + \dfrac{p}{D}\cosh\!\left(\dfrac{R}{\ell}\right)}.
```
The full analytical solution is
```math
    c(x) =
    \frac{\dfrac{p}{D}\,c_{\mathrm{ext}}}
         {\dfrac{1}{\ell}\sinh\!\left(\dfrac{R}{\ell}\right)
          + \dfrac{p}{D}\cosh\!\left(\dfrac{R}{\ell}\right)}
    \cdot \cosh\!\left(\frac{x}{\ell}\right).
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

## Estimating resources

Memory (RAM): usually transferable
Runtime: often very different between computers
Threads scaling: depends on CPU differences
I/O speed: usually very different

### Estimating memory (RAM)
```bash
/usr/bin/time -v python script.py
```

Example for running simulation (with gif and analysis) for simulation with 2500 mesh points.
```bash
User time (seconds): 15.27
System time (seconds): 0.17
Percent of CPU this job got: 98% Elapsed (wall clock) time (h:mm:ss or m:ss): 0:15.67
Average shared text size (kbytes): 0
Average unshared data size (kbytes): 0
Average stack size (kbytes): 0
Average total size (kbytes): 0
Maximum resident set size (kbytes): 192640 # set mem_mb to at least this number / 1024 (for kbytes to mbytes conversion)
Average resident set size (kbytes): 0
Major (requiring I/O) page faults: 0
Minor (reclaiming a frame) page faults: 50329
Voluntary context switches: 3327
Involuntary context switches: 214
Swaps: 0
File system inputs: 0
File system outputs: 1032
Socket messages sent: 0
Socket messages received: 0
Signals delivered: 0
Page size (bytes): 4096
Exit status: 0
```





