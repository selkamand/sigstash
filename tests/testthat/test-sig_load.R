test_that("sig_load works", {
  requireNamespace("sigshared", quietly = TRUE)

  df_available <- sig_available()
  datasets <- df_available[["dataset"]]

  expect_gt(length(datasets), 0)

  # All collections loadable without errors
  for (dataset in datasets) {
    expect_error(sig_load(dataset), NA)
    expect_error(sigshared::assert_signature_collection(sig_load(dataset)), NA)
  }


  # Expect Dataframe Output When format = "tidy"
  expect_s3_class(sig_load(datasets[1], format = "tidy"), "data.frame")

  # Expect Dataframe Output When format = "sigminer"
  expect_s3_class(sig_load(datasets[1], format = "sigminer"), "data.frame")
})


test_that("sig_load_annotations works", {
  requireNamespace("sigshared", quietly = TRUE)

  df_available <- sig_available()
  datasets <- df_available[["dataset"]]

  expect_gt(length(datasets), 0)

  # All collections loadable without errors
  for (dataset in datasets) {
    expect_error(sig_load_annotations(dataset), NA)

    sig <- sig_load(dataset)
    sigs <- sig[["signature"]]
    expect_error(sigshared::assert_signature_annotations(sig_load_annotations(dataset), required_signatures = sigs), NA)
  }
})
