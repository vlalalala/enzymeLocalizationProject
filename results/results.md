# Results

## Optimization
Not a good idea to:
- compare the variance within a round (since there is a strong dependency on the sampling of the trials for that specific round)
- compare the worst vs the best result (since the cutoff condition will then depend more on how well-tuned the search started with)
- 

## Cases without enzymes
### Simple decay $X\rightarrow Y$, maximization of Y
In order to test the optimization procedure, we start with a very simple problem whose analytical solution we know: We study what the optimal location of the semipermeable membrane within the vesicle in order to maximize the flux of species Y. The vesicle is located within a reservoir with a non-zero concentration of species X and zero concentration of species Y. Within the vesicle, there is spontaneous decay of species X into species Y.  

We find that the flux of Y is maximized the closer the membrane is to the origin. This is because there is a concentration jump for X at the membrane, with a larger concentration on the side of the membrane closer to the exterior. Diffusive molecules of X might not pass the semipermeable membrane and potentially escape earlier than in the case that they had not "bumped" against the membrane. These molecules then have less chances to decay into Y.

This begs the question: What happens if Y further decays into Z and we wish to maximize the amoumnt of Y? If we want to maximize the production of Y, the membrane should be as close to the origin as possible. However, the average escape time will increase, making it more likely for Y to decay into Z. The optimal location of the membrane will thus depend on the relation between the decay rates.

### Chained decay $X \rightarrow Y \rightarrow Z$, with one inner boundary
The steady state solution to this problem can be computed analytically. The optimization procedure can thus be checked against the analytical solution.

The idea is to plot the optimal position of the one inner boundary as a function of the decay rates $k_1$ and $k_2$. It makes sense to choose $k_1$ and $k_2$, such that the reaction timescales ($\tau_\mathrm{react} = 1/k$) within ranges spanning from much slower to much faster than the diffusion timescale ($\tau_\mathrm{diff} = R^2 / D$).

If $\tau_\mathrm{diff}/10 < \tau_\mathrm{react} < \tau_\mathrm{diff}*10$, this means
```math
\frac{R^2}{D} \cdot \frac{1}{10} < 1/k < \frac{R^2}{D} \cdot 10 \\
\frac{D}{R^2} \cdot 10 > k > \frac{D}{R^2} \cdot \frac{1}{10}
```
Then, for D=6.6e-11 m^2/s, R=1e-5 m, we get $0.066 1/s < k < 6.6 1/s$ (we do $0.01 1/s < k < 10 1/s$, with logarithmic scales).
