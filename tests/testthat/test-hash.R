test_that("All signatures can be accurately checked using sig_identify_collection", {

  # NOTE: This test is expected to fail after any update or addition to the signature library.
  # If it fails, please run precompute_and_save_md5s() to update the precomputed MD5 sums,
  # and then re-run the tests.

  # Recompute md5sum for all collections
  compute_md5_for_all_collections();

  # Get all available datasets
  available <- sig_available()
  datasets <- available[['dataset']]

  # Define the formats to check
  formats <- c("tidy", "sigstash", "sigminer")

  # Loop over each dataset and format
  for (dataset_name in datasets) {
    for (fmt in formats) {
      # Load the dataset in the current format
      signatures <- sig_load(dataset_name, format = fmt)

      # Check that the result from sig_identify_collection matches the dataset name
      expect_equal(
        sig_identify_collection(signatures, return = "name"),
        dataset_name,
        info = paste("Failed for dataset", dataset_name, "in format", fmt)
      )

      # Check that the result from sig_identify_collection (name + format) matches the expected unique name
      expected_unique_name <- paste(dataset_name, " (format: ", fmt, ")", sep = "")
      expect_equal(
        sig_identify_collection(signatures, return = "name_plus_format"),
        expected_unique_name,
        info = paste("Failed for dataset", dataset_name, "in format", fmt)
      )
    }
  }
})
