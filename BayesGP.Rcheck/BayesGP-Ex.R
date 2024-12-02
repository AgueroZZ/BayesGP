pkgname <- "BayesGP"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "BayesGP-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('BayesGP')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("compute_post_fun_iwp")
### * compute_post_fun_iwp

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: compute_post_fun_iwp
### Title: Computing the posterior samples of the function or its
###   derivative using the posterior samples of the basis coefficients for
###   iwp
### Aliases: compute_post_fun_iwp

### ** Examples

knots <- c(0, 0.2, 0.4, 0.6)
samps <- matrix(rnorm(n = (3 * 10)), ncol = 10)
result <- compute_post_fun_iwp(samps = samps, knots = knots, refined_x = seq(0, 1, by = 0.1), p = 2)
plot(result[, 2] ~ result$x, type = "l", ylim = c(-0.3, 0.3))
for (i in 1:9) {
  lines(result[, (i + 1)] ~ result$x, lty = "dashed", ylim = c(-0.1, 0.1))
}
global_samps <- matrix(rnorm(n = (2 * 10), sd = 0.1), ncol = 10)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("compute_post_fun_iwp", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("compute_weights_precision_helper")
### * compute_weights_precision_helper

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: compute_weights_precision_helper
### Title: Constructing the precision matrix given the knot sequence
###   (helper)
### Aliases: compute_weights_precision_helper

### ** Examples

compute_weights_precision_helper(x = c(0,0.2,0.4,0.6,0.8))



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("compute_weights_precision_helper", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("get_default_option_list_MCMC")
### * get_default_option_list_MCMC

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: get_default_option_list_MCMC
### Title: Get default options for MCMC implementation
### Aliases: get_default_option_list_MCMC

### ** Examples

# Example: Get the default option list
options <- get_default_option_list_MCMC()
print(options)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("get_default_option_list_MCMC", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("global_poly_helper")
### * global_poly_helper

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: global_poly_helper
### Title: Constructing and evaluating the global polynomials, to account
###   for boundary conditions (design matrix)
### Aliases: global_poly_helper

### ** Examples

global_poly_helper(x = c(0, 0.2, 0.4, 0.6, 0.8), p = 2)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("global_poly_helper", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("local_poly_helper")
### * local_poly_helper

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: local_poly_helper
### Title: Constructing and evaluating the local O-spline basis (design
###   matrix)
### Aliases: local_poly_helper

### ** Examples

local_poly_helper(knots = c(0, 0.2, 0.4, 0.6, 0.8), refined_x = seq(0, 0.8, by = 0.1), p = 2)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("local_poly_helper", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
