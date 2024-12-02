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
library(tidyverse)

## -----------------------------------------------------------------------------
head(covid_canada)
par(oldpar)
