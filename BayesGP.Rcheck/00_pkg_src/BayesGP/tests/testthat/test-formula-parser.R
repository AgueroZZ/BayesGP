test_that("multiplication works", {
  expect_equal(2 * 2, 4)
  # Test case 1, only 1 random effect
  parse_result <- parse_formula(rent ~ f(
    smoothing_var = floor.size,
    model = "IWP",
    order = polyOrder1
  ))
  expect_equal(parse_result$response, as.symbol("rent"))
  expect_equal(parse_result$rand_effects[[1]]$smoothing_var, as.symbol("floor.size"))
  expect_equal(parse_result$rand_effects[[1]]$model, "IWP")
  expect_equal(parse_result$rand_effects[[1]]$order, as.symbol("polyOrder1"))

  # Test case 2, 2 random effects
  parse_result <- parse_formula(rent ~ f(
    smoothing_var = floor.size,
    model = "IWP",
    order = polyOrder1
  )
  + f(
      smoothing_var = year,
      model = "IWP",
      order = polyOrder2, k = 10, # should add a checker for k >= 3
      sd.prior = list(prior = "exp", para = list(u = 1, alpha = 0.5)),
      boundary.prior = list(prec = 0.01)
    ))
  expect_equal(parse_result$response, as.symbol("rent"))
  expect_equal(parse_result$rand_effects[[1]]$smoothing_var, as.symbol("floor.size"))
  expect_equal(parse_result$rand_effects[[1]]$model, "IWP")
  expect_equal(parse_result$rand_effects[[1]]$order, as.symbol("polyOrder1"))
  expect_equal(parse_result$rand_effects[[2]]$smoothing_var, as.symbol("year"))
  expect_equal(parse_result$rand_effects[[2]]$model, "IWP")
  expect_equal(parse_result$rand_effects[[2]]$order, as.symbol("polyOrder2"))
  expect_equal(parse_result$rand_effects[[2]]$k, 10)

  # Test case 3, 2 random effects, 2 fixed effects
  parse_result <- parse_formula(rent ~ location + f(
    smoothing_var = floor.size,
    model = "IWP",
    order = polyOrder1
  )
  + score + f(
      smoothing_var = year,
      model = "IWP",
      order = polyOrder2, k = 10, # should add a checker for k >= 3
      sd.prior = list(prior = "exp", para = list(u = 1, alpha = 0.5)),
      boundary.prior = list(prec = 0.01)
    ))
  expect_equal(parse_result$response, as.symbol("rent"))
  expect_equal(parse_result$rand_effects[[1]]$smoothing_var, as.symbol("floor.size"))
  expect_equal(parse_result$rand_effects[[1]]$model, "IWP")
  expect_equal(parse_result$rand_effects[[1]]$order, as.symbol("polyOrder1"))
  expect_equal(parse_result$rand_effects[[2]]$smoothing_var, as.symbol("year"))
  expect_equal(parse_result$rand_effects[[2]]$model, "IWP")
  expect_equal(parse_result$rand_effects[[2]]$order, as.symbol("polyOrder2"))
  expect_equal(parse_result$rand_effects[[2]]$k, 10)
  expect_equal(parse_result$fixed_effects[[1]], as.symbol("location"))
  expect_equal(parse_result$fixed_effects[[2]], as.symbol("score"))
})
