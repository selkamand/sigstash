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


cosmic_style_indels <- c("1:Del:C:0", "1:Del:C:1", "1:Del:C:2", "1:Del:C:3", "1:Del:C:4",
                         "1:Del:C:5", "1:Del:T:0", "1:Del:T:1", "1:Del:T:2", "1:Del:T:3",
                         "1:Del:T:4", "1:Del:T:5", "1:Ins:C:0", "1:Ins:C:1", "1:Ins:C:2",
                         "1:Ins:C:3", "1:Ins:C:4", "1:Ins:C:5", "1:Ins:T:0", "1:Ins:T:1",
                         "1:Ins:T:2", "1:Ins:T:3", "1:Ins:T:4", "1:Ins:T:5", "2:Del:R:0",
                         "2:Del:R:1", "2:Del:R:2", "2:Del:R:3", "2:Del:R:4", "2:Del:R:5",
                         "3:Del:R:0", "3:Del:R:1", "3:Del:R:2", "3:Del:R:3", "3:Del:R:4",
                         "3:Del:R:5", "4:Del:R:0", "4:Del:R:1", "4:Del:R:2", "4:Del:R:3",
                         "4:Del:R:4", "4:Del:R:5", "5:Del:R:0", "5:Del:R:1", "5:Del:R:2",
                         "5:Del:R:3", "5:Del:R:4", "5:Del:R:5", "2:Ins:R:0", "2:Ins:R:1",
                         "2:Ins:R:2", "2:Ins:R:3", "2:Ins:R:4", "2:Ins:R:5", "3:Ins:R:0",
                         "3:Ins:R:1", "3:Ins:R:2", "3:Ins:R:3", "3:Ins:R:4", "3:Ins:R:5",
                         "4:Ins:R:0", "4:Ins:R:1", "4:Ins:R:2", "4:Ins:R:3", "4:Ins:R:4",
                         "4:Ins:R:5", "5:Ins:R:0", "5:Ins:R:1", "5:Ins:R:2", "5:Ins:R:3",
                         "5:Ins:R:4", "5:Ins:R:5", "2:Del:M:1", "3:Del:M:1", "3:Del:M:2",
                         "4:Del:M:1", "4:Del:M:2", "4:Del:M:3", "5:Del:M:1", "5:Del:M:2",
                         "5:Del:M:3", "5:Del:M:4", "5:Del:M:5")

sigminer_style_indels <- c("1_Del_C_0", "1_Del_C_1", "1_Del_C_2", "1_Del_C_3", "1_Del_C_4",
                          "1_Del_C_5", "1_Del_T_0", "1_Del_T_1", "1_Del_T_2", "1_Del_T_3",
                          "1_Del_T_4", "1_Del_T_5", "1_Ins_C_0", "1_Ins_C_1", "1_Ins_C_2",
                          "1_Ins_C_3", "1_Ins_C_4", "1_Ins_C_5", "1_Ins_T_0", "1_Ins_T_1",
                          "1_Ins_T_2", "1_Ins_T_3", "1_Ins_T_4", "1_Ins_T_5", "2_Del_R_0",
                          "2_Del_R_1", "2_Del_R_2", "2_Del_R_3", "2_Del_R_4", "2_Del_R_5",
                          "3_Del_R_0", "3_Del_R_1", "3_Del_R_2", "3_Del_R_3", "3_Del_R_4",
                          "3_Del_R_5", "4_Del_R_0", "4_Del_R_1", "4_Del_R_2", "4_Del_R_3",
                          "4_Del_R_4", "4_Del_R_5", "5_Del_R_0", "5_Del_R_1", "5_Del_R_2",
                          "5_Del_R_3", "5_Del_R_4", "5_Del_R_5", "2_Ins_R_0", "2_Ins_R_1",
                          "2_Ins_R_2", "2_Ins_R_3", "2_Ins_R_4", "2_Ins_R_5", "3_Ins_R_0",
                          "3_Ins_R_1", "3_Ins_R_2", "3_Ins_R_3", "3_Ins_R_4", "3_Ins_R_5",
                          "4_Ins_R_0", "4_Ins_R_1", "4_Ins_R_2", "4_Ins_R_3", "4_Ins_R_4",
                          "4_Ins_R_5", "5_Ins_R_0", "5_Ins_R_1", "5_Ins_R_2", "5_Ins_R_3",
                          "5_Ins_R_4", "5_Ins_R_5", "2_Del_M_1", "3_Del_M_1", "3_Del_M_2",
                          "4_Del_M_1", "4_Del_M_2", "4_Del_M_3", "5_Del_M_1", "5_Del_M_2",
                          "5_Del_M_3", "5_Del_M_4", "5_Del_M_5")

cosmic_style_sv <- c("clustered_del_1-10Kb", "clustered_del_10-100Kb", "clustered_del_100Kb-1Mb",
                     "clustered_del_1Mb-10Mb", "clustered_del_>10Mb", "clustered_tds_1-10Kb",
                     "clustered_tds_10-100Kb", "clustered_tds_100Kb-1Mb", "clustered_tds_1Mb-10Mb",
                     "clustered_tds_>10Mb", "clustered_inv_1-10Kb", "clustered_inv_10-100Kb",
                     "clustered_inv_100Kb-1Mb", "clustered_inv_1Mb-10Mb", "clustered_inv_>10Mb",
                     "clustered_trans", "non-clustered_del_1-10Kb", "non-clustered_del_10-100Kb",
                     "non-clustered_del_100Kb-1Mb", "non-clustered_del_1Mb-10Mb",
                     "non-clustered_del_>10Mb", "non-clustered_tds_1-10Kb", "non-clustered_tds_10-100Kb",
                     "non-clustered_tds_100Kb-1Mb", "non-clustered_tds_1Mb-10Mb",
                     "non-clustered_tds_>10Mb", "non-clustered_inv_1-10Kb", "non-clustered_inv_10-100Kb",
                     "non-clustered_inv_100Kb-1Mb", "non-clustered_inv_1Mb-10Mb",
                     "non-clustered_inv_>10Mb", "non-clustered_trans")

sigminer_style_sv <- c("clustered:del:1-10Kb", "clustered:del:10-100Kb", "clustered:del:100Kb-1Mb",
                       "clustered:del:1Mb-10Mb", "clustered:del:>10Mb", "clustered:tds:1-10Kb",
                       "clustered:tds:10-100Kb", "clustered:tds:100Kb-1Mb", "clustered:tds:1Mb-10Mb",
                       "clustered:tds:>10Mb", "clustered:inv:1-10Kb", "clustered:inv:10-100Kb",
                       "clustered:inv:100Kb-1Mb", "clustered:inv:1Mb-10Mb", "clustered:inv:>10Mb",
                       "clustered:trans", "non-clustered:del:1-10Kb", "non-clustered:del:10-100Kb",
                       "non-clustered:del:100Kb-1Mb", "non-clustered:del:1Mb-10Mb",
                       "non-clustered:del:>10Mb", "non-clustered:tds:1-10Kb", "non-clustered:tds:10-100Kb",
                       "non-clustered:tds:100Kb-1Mb", "non-clustered:tds:1Mb-10Mb",
                       "non-clustered:tds:>10Mb", "non-clustered:inv:1-10Kb", "non-clustered:inv:10-100Kb",
                       "non-clustered:inv:100Kb-1Mb", "non-clustered:inv:1Mb-10Mb",
                       "non-clustered:inv:>10Mb", "non-clustered:trans")

cosmic_style_cn48 <- sig_get_valid_cosmic_channels("CN48")
sigminer_style_cn48 <- c("0_homdel_0-100Kb", "0_homdel_100Kb-1Mb", "0_homdel_>1Mb",
                         "1_LOH_0-100Kb", "1_LOH_100Kb-1Mb", "1_LOH_1Mb-10Mb", "1_LOH_10Mb-40Mb",
                         "1_LOH_>40Mb", "2_LOH_0-100Kb", "2_LOH_100Kb-1Mb", "2_LOH_1Mb-10Mb",
                         "2_LOH_10Mb-40Mb", "2_LOH_>40Mb", "3-4_LOH_0-100Kb", "3-4_LOH_100Kb-1Mb",
                         "3-4_LOH_1Mb-10Mb", "3-4_LOH_10Mb-40Mb", "3-4_LOH_>40Mb", "5-8_LOH_0-100Kb",
                         "5-8_LOH_100Kb-1Mb", "5-8_LOH_1Mb-10Mb", "5-8_LOH_10Mb-40Mb",
                         "5-8_LOH_>40Mb", "9+_LOH_0-100Kb", "9+_LOH_100Kb-1Mb", "9+_LOH_1Mb-10Mb",
                         "9+_LOH_10Mb-40Mb", "9+_LOH_>40Mb", "2_het_0-100Kb", "2_het_100Kb-1Mb",
                         "2_het_1Mb-10Mb", "2_het_10Mb-40Mb", "2_het_>40Mb", "3-4_het_0-100Kb",
                         "3-4_het_100Kb-1Mb", "3-4_het_1Mb-10Mb", "3-4_het_10Mb-40Mb",
                         "3-4_het_>40Mb", "5-8_het_0-100Kb", "5-8_het_100Kb-1Mb", "5-8_het_1Mb-10Mb",
                         "5-8_het_10Mb-40Mb", "5-8_het_>40Mb", "9+_het_0-100Kb", "9+_het_100Kb-1Mb",
                         "9+_het_1Mb-10Mb", "9+_het_10Mb-40Mb", "9+_het_>40Mb")

cosmic_style_sbs96 <- sig_get_valid_cosmic_channels("SBS96")
sigminer_style_sbs96 <- c("A[C>A]A", "A[C>A]C", "A[C>A]G", "A[C>A]T", "A[C>G]A", "A[C>G]C",
                          "A[C>G]G", "A[C>G]T", "A[C>T]A", "A[C>T]C", "A[C>T]G", "A[C>T]T",
                          "A[T>A]A", "A[T>A]C", "A[T>A]G", "A[T>A]T", "A[T>C]A", "A[T>C]C",
                          "A[T>C]G", "A[T>C]T", "A[T>G]A", "A[T>G]C", "A[T>G]G", "A[T>G]T",
                          "C[C>A]A", "C[C>A]C", "C[C>A]G", "C[C>A]T", "C[C>G]A", "C[C>G]C",
                          "C[C>G]G", "C[C>G]T", "C[C>T]A", "C[C>T]C", "C[C>T]G", "C[C>T]T",
                          "C[T>A]A", "C[T>A]C", "C[T>A]G", "C[T>A]T", "C[T>C]A", "C[T>C]C",
                          "C[T>C]G", "C[T>C]T", "C[T>G]A", "C[T>G]C", "C[T>G]G", "C[T>G]T",
                          "G[C>A]A", "G[C>A]C", "G[C>A]G", "G[C>A]T", "G[C>G]A", "G[C>G]C",
                          "G[C>G]G", "G[C>G]T", "G[C>T]A", "G[C>T]C", "G[C>T]G", "G[C>T]T",
                          "G[T>A]A", "G[T>A]C", "G[T>A]G", "G[T>A]T", "G[T>C]A", "G[T>C]C",
                          "G[T>C]G", "G[T>C]T", "G[T>G]A", "G[T>G]C", "G[T>G]G", "G[T>G]T",
                          "T[C>A]A", "T[C>A]C", "T[C>A]G", "T[C>A]T", "T[C>G]A", "T[C>G]C",
                          "T[C>G]G", "T[C>G]T", "T[C>T]A", "T[C>T]C", "T[C>T]G", "T[C>T]T",
                          "T[T>A]A", "T[T>A]C", "T[T>A]G", "T[T>A]T", "T[T>C]A", "T[T>C]C",
                          "T[T>C]G", "T[T>C]T", "T[T>G]A", "T[T>G]C", "T[T>G]G", "T[T>G]T"
)
cosmic_style_dbs78 <- sig_get_valid_cosmic_channels("DBS78")
sigminer_style_dbs78 <- c("AC>CA", "AC>CG", "AC>CT", "AC>GA", "AC>GG", "AC>GT", "AC>TA",
                          "AC>TG", "AC>TT", "AT>CA", "AT>CC", "AT>CG", "AT>GA", "AT>GC",
                          "AT>TA", "CC>AA", "CC>AG", "CC>AT", "CC>GA", "CC>GG", "CC>GT",
                          "CC>TA", "CC>TG", "CC>TT", "CG>AT", "CG>GC", "CG>GT", "CG>TA",
                          "CG>TC", "CG>TT", "CT>AA", "CT>AC", "CT>AG", "CT>GA", "CT>GC",
                          "CT>GG", "CT>TA", "CT>TC", "CT>TG", "GC>AA", "GC>AG", "GC>AT",
                          "GC>CA", "GC>CG", "GC>TA", "TA>AT", "TA>CG", "TA>CT", "TA>GC",
                          "TA>GG", "TA>GT", "TC>AA", "TC>AG", "TC>AT", "TC>CA", "TC>CG",
                          "TC>CT", "TC>GA", "TC>GG", "TC>GT", "TG>AA", "TG>AC", "TG>AT",
                          "TG>CA", "TG>CC", "TG>CT", "TG>GA", "TG>GC", "TG>GT", "TT>AA",
                          "TT>AC", "TT>AG", "TT>CA", "TT>CC", "TT>CG", "TT>GA", "TT>GC",
                          "TT>GG")

test_that("Channel to Type conversion works", {
  expect_no_error(sig_convert_channel2type(sbs1536_channels, sigclass = "SBS1536"))
  expect_equal(unname(sig_convert_channel2type(sbs1536_channels, sigclass = "SBS1536")), sbs1536_type)

  expect_error(sig_convert_channel2type(sbs1536_channels, sigclass = "SBS96"))

  # Sigclass Chass
  expect_no_error(sig_valid_sigclass())

  for (sigclass in sig_valid_sigclass()) {
    expect_no_error(sig_get_valid_cosmic_channels(sigclass))
    expect_no_error(sig_get_valid_cosmic_types(sigclass))
  }
})


test_that("Cosmic to Sigminer conversion works", {
    expect_no_error(sig_convert_channel_name(cosmic_style_indels, from = "cosmic", to="sigminer"))

  # Indels
  expect_equal(
      sig_convert_channel_name(cosmic_style_indels, from = "cosmic", to="sigminer"),
      sigminer_style_indels
    )
    expect_equal(
      sig_convert_channel_name(sigminer_style_indels, from = "sigminer", to="cosmic"),
      cosmic_style_indels
    )

    # SVs
    expect_equal(
      sig_convert_channel_name(cosmic_style_sv, from = "cosmic", to="sigminer"),
      sigminer_style_sv
    )
    expect_equal(
      sig_convert_channel_name(sigminer_style_sv, from = "sigminer", to="cosmic"),
      cosmic_style_sv
    )

    # Copynumber
    expect_equal(
      sig_convert_channel_name(cosmic_style_cn48, from = "cosmic", to="sigminer"),
      sigminer_style_cn48
    )
    expect_equal(
      sig_convert_channel_name(sigminer_style_cn48, from = "sigminer", to="cosmic"),
      cosmic_style_cn48
    )

    # SBS96
    expect_equal(
      sig_convert_channel_name(cosmic_style_sbs96, from = "cosmic", to="sigminer"),
      sigminer_style_sbs96
    )
    expect_equal(
      sig_convert_channel_name(sigminer_style_sbs96, from = "sigminer", to="cosmic"),
      cosmic_style_sbs96
    )

    # DBS78
    expect_equal(
      sig_convert_channel_name(cosmic_style_dbs78, from = "cosmic", to="sigminer"),
      sigminer_style_dbs78
    )
    expect_equal(
      sig_convert_channel_name(sigminer_style_dbs78, from = "sigminer", to="cosmic"),
      cosmic_style_dbs78
    )
})
