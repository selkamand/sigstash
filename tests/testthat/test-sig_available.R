test_that("sig_available", {

  # Runs without error
  expect_error(sig_available(), NA)


  df_available = sig_available()

  # Class is data.frame
  expect_s3_class(df_available, "data.frame")

  # Expected colnames
  expect_equal(c('dataset', 'description'), colnames(df_available), ignore_attr = TRUE)

  # Data.frame has at least 1 row
  expect_true(nrow(df_available) > 0)
})
