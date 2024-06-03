path_inst <- system.file(package = "sigstash")

# vec_file_sigclass <- c(
#   "COSMIC_v3.4_CN_GRCh37.txt" = "CN",
#   "COSMIC_v3.4_DBS_GRCh37.txt" = "DBS",
#   "COSMIC_v3.4_DBS_GRCh38.txt" = "DBS",
#   "COSMIC_v3.4_SBS_GRCh37.txt" = "SBS",
#   "COSMIC_v3.4_SBS_GRCh38.txt" = "SBS",
#   "COSMIC_v3.4_ID_GRCh37.txt" = "ID",
#   "COSMIC_v3.4_SV_GRCh38.txt" = "SV",
#   "COSMIC_v3.4_RNA-SBS_GRCh37.txt" = "RNA-SBS"
# )


for (i in seq_along(vec_file_sigclass)) {
  orig_file <- names(vec_file_sigclass)[i]
  sigclass <- vec_file_sigclass[i]

  path <- system.file("original_files/cosmic", orig_file, package = "sigstash")
  data <- utils::read.csv(path, header = TRUE, sep = "\t")

  # Convert to sigstash-style
  df_sigs <- sig_cosmic_to_sigstash(data, sigclass = sigclass)

  # Write output
  outfile <- paste0(tools::file_path_sans_ext(basename(orig_file)), ".csv")
  message("Writing: ", outfile)

  utils::write.csv(
    df_sigs,
    file = paste0(path_inst, "/reference_signatures/", outfile),
    row.names = FALSE,
    quote = TRUE
  )
}
