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
