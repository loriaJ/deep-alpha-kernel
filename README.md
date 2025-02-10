# Deep Kernel Posterior Learning under Infinite Variance Prior Weights 

**"Deep Kernel Posterior Learning under Infinite Variance Prior Weights" Authors: [Jorge Loría](https://loriaj.github.io/), [Anindya Bhadra](https://www.stat.purdue.edu/~bhadra/).** 

Preprint: [https://arxiv.org/abs/2410.01284](https://arxiv.org/abs/2410.01284)

To replicate the simulations use the files: 
- `example_1d.R`, and
- `example_2d.R`.
The library `mvtnorm` is required to run the main files. The example for two dimensions needs `ggplot2` for visualizing the results.

The main functions to run the predictions are in: `function_several_layers.R` file. Specifically, the main function used is: `stable_KP`. This function takes 3 required parameters:

- `x0`: the observed locations matrix,
- `y0`: the observed responses, and
- `x_new`: the locations for predictions.
 
It can also receive the optional parameters:
- `n_iters`, the number of simulations, predefined to be 3000, must be a positive integer. 
- `alpha2`, the index parameter of the $\alpha$-stables, predefined to be 1, must be in the open interval $(0,2)$.
- `sigma_offset`, the assumed standard deviation, predefined to be 1, must be positive. Only used if `sigma_known=TRUE`.
- `sigma_known` if the assumed standard deviation is given by the `sigma_offset` parameter, or `FALSE` if it is not known, predefined to be `FALSE`,
- `n_layers`, number of hidden layers, default value is 2.
- `sd_offset`, scale of the half-Cauchy prior of the standard deviation.
- `anisotropic`, if TRUE (default) uses the $J_\delta(\theta)$ established in Th. 1 of the main paper. If FALSE, need to specify a function $k_0$, in the outer environment that will act as a kernel function, e.g. Matérn, RBF, etc.  
- `delta`, power of the activation function $g_\delta(\xi)=\xi^\delta \mathbf{1}_{\{\xi >0\}}$, must be a non-negative integer. Default value of 1.
- `save_sigmas_indic`, if the posterior covariance matrices should be saved or noted. Default is not to save. When saving use caution as it requires more RAM.


