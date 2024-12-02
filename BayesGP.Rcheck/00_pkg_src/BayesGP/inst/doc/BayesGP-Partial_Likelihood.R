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
data <- as.data.frame(ccData)
data$exposure <- data$exposure
mod <- model_fit(formula = case ~ f(x = exposure, 
                                    model = "IWP", 
                                    order = 2, k = 30,
                                    initial_location = median(data$exposure), 
                                    sd.prior = list(prior = "exp", param = list(u = 1, alpha = 0.5), h = 1)),
                 family = "cc",
                 strata = "subject",
                 weight = NULL,
                 data = data,
                 method = "aghq")


## -----------------------------------------------------------------------------
true_effect <- function(x) {3 *(x^2 - .5^2)}
plot(mod)
lines(I(true_effect(seq(0,1,by = 0.1)) - true_effect(median(data$exposure))) ~ seq(0,1,by = 0.1), col = "red")

## -----------------------------------------------------------------------------
data <- survival::kidney
head(data)
mod <- model_fit(formula = time ~ age + sex + f(x = id, 
                                    model = "IID", 
                                    sd.prior = list(prior = "exp", param = list(u = 1, alpha = 0.5))),
                 family = "coxph",
                 cens = "status",
                 data = data,
                 method = "aghq")

## -----------------------------------------------------------------------------
samps_age <- sample_fixed_effect(mod, variables = "age")
samps_sex <- sample_fixed_effect(mod, variables = "sex")
par(mfrow = c(1,2))
hist(samps_age, main = "Samples for effect of age", xlab = "Effect")
hist(samps_sex, main = "Samples for effect of sex", xlab = "Effect")
par(oldpar)
