test_that("surv.comp function works", {

  # load dataset
  data("EarlyEffectData")

  # run your function
  result <- surv.comp(
    time = EarlyEffectData$time,
    status = EarlyEffectData$status,
    group = EarlyEffectData$group,
    n_perm = 10,
    censor_type = "right",
    method = "joint_method"
  )

  # Check class
  expect_s3_class(result, "survival_comp_roc")

  # Check list elements
  expect_true("message" %in% names(result))
  expect_true("result" %in% names(result))

  # Check message is character
  expect_type(result$message, "character")

  # Check result is data frame
  expect_s3_class(result$result, "data.frame")

  # Check columns exist
  expect_true(all(c("estimate", "p_value") %in% colnames(result$result)))

  # Check p-value range
  expect_true(all(result$result$p_value >= 0))
  expect_true(all(result$result$p_value <= 1))

})
