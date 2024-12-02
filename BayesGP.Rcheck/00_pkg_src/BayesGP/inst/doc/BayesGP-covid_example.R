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
head(covid_canada)

## ----warning=FALSE------------------------------------------------------------
fit_result <- model_fit(new_deaths ~ weekdays1 + weekdays2 + weekdays3 + weekdays4 + weekdays5 + weekdays6 +
                          f(smoothing_var = t, model = "IWP", order = 3, k = 100, sd.prior = list(prior = "exp", param = list(u = 0.02, alpha = 0.5), h = 1)), 
                        data = covid_canada, method = "aghq", family = "Poisson")

## -----------------------------------------------------------------------------
summary(fit_result)

## ----warning=FALSE------------------------------------------------------------
plot(fit_result)

## ----warning=FALSE------------------------------------------------------------
predict_f <- predict(fit_result, variable = "t", newdata = data.frame(t = seq(from = 605, to = 615, by = 0.1)))
matplot(x = predict_f[,1], y = predict_f[,2:4], type = 'l', ylab = "f", xlab = "t", col = c("grey", "grey", "black"), lty = c("dashed", "dashed", "solid"))

## ----warning=FALSE------------------------------------------------------------
predict_f1st <- predict(fit_result, variable = "t", newdata = data.frame(t = seq(from = 605, to = 615, by = 0.1)), deriv = 1)
matplot(x = predict_f1st[,1], y = predict_f1st[,2:4], type = 'l', ylab = "f'", xlab = "t", col = c("grey", "grey", "black"), lty = c("dashed", "dashed", "solid"))

## ----warning=FALSE------------------------------------------------------------
predict_f2nd <- predict(fit_result, variable = "t", newdata = data.frame(t = seq(from = 605, to = 617, by = 0.1)), deriv = 2)
matplot(x = predict_f2nd[,1], y = predict_f2nd[,2:4], type = 'l', ylab = "f''", xlab = "t", col = c("grey", "grey", "black"), lty = c("dashed", "dashed", "solid"))
par(oldpar)
