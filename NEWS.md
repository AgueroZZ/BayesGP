# BayesGP 0.1.2 (2024-06-16)

* Bug fix: 
- Fixed a bug in the `global_poly_helper_sgp` function and subsequently in the `predict` method that caused the prediction to fail when the smoothing variable is not pre-initialized manually. This should not affect models and predictions with manually initialized smoothing variables.

* Minor changes:
- The definitions of `initial_location` being `left`, `right` and `middle` in IWP and sGP are not defined based on the region rather than the provided smoothing variable.

# BayesGP 0.1.3 (2024-10-01)

* Bug fix:
- Fixed a bug in the structure of posterior samples when the model_fit object is of class "nlmb", now the object has posterior samples stored as a list `object$samps$samps`, which prevents the summary method from unexpected crashes.

* Minor changes:
- `var_plot` and `var_density` are now called `sd_plot` and `sd_density`.
- In model_fit, `knots` can now take a single numeric value, which will then be assumed as `k`.
- The default `order` of iwp is assumed to be 2.
- Plot now produce multiple figures of both the smoothing results and the SD parameters, which could be proceeded by pressing "next".


# BayesGP 0.1.4 (2024-12-03)

* Minor changes:
- The initial location in `f` now could also support input of `mean`, `min` and `max`.
- Update the function `sample_fixed_effect` to accept an input of `M`, denoting the number of posterior samples to take.
- Update the method `plot`, to accept an argument of `select`, which selects which variable to show in the plot.
- Create functions `sample_sd` that is used to extract posterior samples of the sd parameter, as well as its subsequent wrapper function `para_samps` which extracts posterior samples of any parameter.
- `custom_template` now takes input in lower case rather than in upper case.




