test_that("sig_write works", {
  signatures <- list(
    sig1 = tibble::tibble(
      type = rep("A>G", 3L),
      channel = c("A[A->G]G", "A[A->G]C", "A[A->G]T"),
      fraction = c(0.4, 0.1, 0.5)
    ),
    sig2 = tibble::tibble(
      type = rep("A>G", 3L),
      channel = c("A[A->G]G", "A[A->G]C", "A[A->G]T"),
      fraction = c(0.4, 0.1, 0.5)
    )
  )

  # RDS
  outfile_rds = withr::local_file(.file = "test.Rds")
  sig_write_signatures(signatures = signatures, filepath = outfile_rds, format = "rds")
  expect_equal(sig_read_signatures(outfile_rds, format = "rds"), signatures)

  # Tidy csv
  outfile_tidy = withr::local_file(.file = "test.tidy.csv")
  sig_write_signatures(signatures = signatures, filepath = outfile_tidy, format = "csv_tidy")
  sig_parsed <- expect_equal(sig_read_signatures(filepath = outfile_tidy, format = "csv_tidy"), signatures, ignore_attr = TRUE)

  # Test attribute addition works as expected
  expect_equal(attr(sig_parsed, "format"), "sigstash")
  expect_equal(attr(sig_parsed, "collection_name"), "test")

  # Test setting custom collection name works as expected
  sig_parsed2 <- sig_read_signatures(filepath = outfile_tidy, format = "csv_tidy", collection_name = "mycollection")
  expect_equal(attr(sig_parsed2, "format"), "sigstash")
  expect_equal(attr(sig_parsed2, "collection_name"), "mycollection")


  # Wide CSV
  outfile_wide = withr::local_file(.file = "test.wide.csv")
  sig_write_signatures(signatures = signatures, filepath = outfile_wide, format = "csv_wide")
  expect_snapshot_file(outfile_wide)
})
