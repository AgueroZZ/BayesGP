oldpar <- par(no.readonly = TRUE)  # Save current settings
## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>", fig.height=3, fig.width=5, margins=TRUE
)
knitr::knit_hooks$set(margins = function(before, options, envir) {
  if (!before) return()
  graphics::par(mar = c(1.5 + 0.9, 1.5 + 0.9, 0.2, 0.2), mgp = c(1.45, 0.45, 0), cex = 1.25, bty='n')
})

## ----setup--------------------------------------------------------------------

library(BayesGP)

## -----------------------------------------------------------------------------
data <- data.frame(year = seq(1821, 1934, by = 1), y = as.numeric(lynx))
data$x <- data$year - min(data$year)
plot(lynx)

## -----------------------------------------------------------------------------
prior_PSD <- list(u = 1, alpha = 0.01)
prior_SD <- BayesGP::prior_conversion_sgp(d = 50, prior = prior_PSD, a = 2*pi/10)

## -----------------------------------------------------------------------------
mod <- model_fit(formula = y ~ f(x = year, 
                                    model = "sGP", 
                                    a = 2*pi/10, k = 30,
                                    sd.prior = list(prior = "exp", param = prior_SD, h = 2)) + 
                               f(x = x, 
                                    model = "IID", 
                                    sd.prior = list(prior = "exp", param = 0.5)),
                   
                 family = "Poisson",
                 data = data,
                 method = "aghq")

## -----------------------------------------------------------------------------
summary(mod)

## -----------------------------------------------------------------------------
pred_sGP <- predict(mod, variable = "year", newdata = data.frame(year = seq(min(data$year), max(data$year), by = 0.1)))
plot(mean ~ year, data = pred_sGP, type = 'l', col = 'black')
lines(q0.025 ~ year, data = pred_sGP, lty = 'dashed', col = 'grey')
lines(q0.975 ~ year, data = pred_sGP, lty = 'dashed', col = 'grey')

## -----------------------------------------------------------------------------
mod <- model_fit(formula = y ~ f(x = year, 
                                    model = "sGP", 
                                    a = 2*pi/10, k = 30,
                                    sd.prior = list(prior = "exp", param = prior_SD, h = 2),
                                    boundary.prior = list(prec = c(Inf, Inf))) + 
                               f(x = x, 
                                    model = "IID", 
                                    sd.prior = list(prior = "exp", param = 0.5)),
                   
                 family = "Poisson",
                 data = data,
                 method = "aghq")
par(oldpar)
