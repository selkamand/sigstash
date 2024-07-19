sbs1536_channels <- c(
  "AA[C>A]GC", "TA[C>G]TA", "TT[C>G]TG", "AC[T>G]GG", "TA[C>G]CG",
  "TT[T>C]CC", "TA[C>A]GT", "AA[T>G]TA", "CC[C>G]TT", "CG[T>G]GG",
  "TA[T>G]CA", "AG[C>G]AA", "AA[C>T]TA", "CA[C>A]GC", "TA[C>T]GC",
  "TC[C>G]TG", "GG[C>G]AC", "TT[T>G]TT", "AT[T>C]TA", "GA[T>G]AA",
  "CC[T>A]TG", "AC[C>T]AA", "GA[T>A]AT", "GC[T>C]TG", "TT[C>T]CG",
  "CC[T>C]AA", "AA[C>G]GT", "GC[C>T]TG", "CC[C>A]GA", "AT[C>T]GA",
  "GC[T>A]GT", "TG[T>A]GG", "TT[T>G]TA", "TG[T>C]TG", "CA[C>T]AC",
  "GT[C>G]TT", "AC[C>G]CA", "GA[C>G]TT", "TT[T>G]CT", "GC[T>C]AT",
  "GG[C>A]GG", "GC[T>G]GC", "TT[C>T]GC", "AC[T>G]AA", "CA[T>A]TA",
  "AA[C>G]AG", "CG[C>G]AG", "AG[T>C]GA", "AT[C>A]TG", "TA[T>G]GC"
)
sbs1536_type <- c(
  "C>A", "C>G", "C>G", "T>G", "C>G", "T>C", "C>A", "T>G", "C>G",
  "T>G", "T>G", "C>G", "C>T", "C>A", "C>T", "C>G", "C>G", "T>G",
  "T>C", "T>G", "T>A", "C>T", "T>A", "T>C", "C>T", "T>C", "C>G",
  "C>T", "C>A", "C>T", "T>A", "T>A", "T>G", "T>C", "C>T", "C>G",
  "C>G", "C>G", "T>G", "T>C", "C>A", "T>G", "C>T", "T>G", "T>A",
  "C>G", "C>G", "T>C", "C>A", "T>G"
)


test_that("Channel to Type conversion works", {
  expect_no_error(sig_convert_channel2type(sbs1536_channels, sigclass = "SBS1536"))
  expect_equal(unname(sig_convert_channel2type(sbs1536_channels, sigclass = "SBS1536")), sbs1536_type)

  expect_error(sig_convert_channel2type(sbs1536_channels, sigclass = "SBS96"))

  # Sigclass Chass
  expect_no_error(sig_valid_sigclass())

  for (sigclass in sig_valid_sigclass()){
    expect_no_error(sig_get_valid_cosmic_channels(sigclass))
    expect_no_error(sig_get_valid_cosmic_types(sigclass))
  }
})
