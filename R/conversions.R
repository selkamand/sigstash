# Convert: COSMIC > SIGVERSE DataFrame
#
# Cosmic Style to sigverse style signature collection files.
# For example of use see inst/original_files/cosmic/comsic_to_sigverse.R script
#
# @param data cosmic style dataframe.
# Has header line, 1 column for each signature, and first column corresponds to the channel.
# Each Sig column is a numeric vector describing the fraction of mutations in each group
#
# @return list of data.frames
#
# @examples
#
# orig_file = 'COSMIC_v3.3.1_SBS_GRCh37.txt'
# path = system.file("original_files/cosmic", orig_file, package = "sigstash")
# data = utils::read.csv(path, header = TRUE, sep = "\t")
#
# # Convert to sigstash-style
# df_sigs = sig_cosmic_to_sigstash(data, sigclass = 'SNV')
#
# # Write output
# outfile = paste0(tools::file_path_sans_ext(basename(orig_file)), '.csv')
# utils::write.csv(df_sigs, file = paste0("inst/reference_signatures/", outfile), row.names = FALSE, quote = TRUE)
sig_cosmic_to_sigstash <- function(data, sigclass = c("SBS", "ID", "CN", "DBS", "SV", "RNA-SBS")) {
  requireNamespace("rlang", quietly = TRUE)

  sigclass <- rlang::arg_match(sigclass)

  signames <- colnames(data)[-1]
  channels <- data[[1]]
  types <- channel2type(channels, sigclass)

  # For each signature ...
  ls_signatures <- lapply(X = signames, FUN = \(sig){
    # Create a data.frame with just Type (renamed to channel) and fractions
    df_sig <- data[, c("Type", sig)]
    colnames(df_sig) <- c("channel", "fraction")

    # Add sigverse defined 'type' column, and 'signature' column describing name
    df_sig[["type"]] <- types
    df_sig[["signature"]] <- sig
    df_sig <- df_sig[, c("signature", "type", "channel", "fraction")]

    # Convert to tibble
    df_sig <- tibble::tibble(df_sig)

    # Return 3 column tibble (channel, type, fraction)
    return(df_sig)
  })

  df_signatures <- do.call("rbind", ls_signatures)

  return(df_signatures)
  # names(ls_signatures) <- signames

  # return(ls_signatures)
}

channel2type <- function(channel, sigclass = c("SBS", "ID", "CN", "DBS", "SV", "RNA-SBS")) {
  requireNamespace("rlang", quietly = TRUE)
  sigclass <- rlang::arg_match(sigclass)

  if (sigclass == "SBS") {
    vec_channel2type <- cosmic_sbs_channel_to_type()
  } else if (sigclass == "CN") {
    vec_channel2type <- cosmic_cn_channel_to_type()
  } else if (sigclass == "ID") {
    vec_channel2type <- cosmic_id_channel_to_type()
  } else if (sigclass == "DBS") {
    vec_channel2type <- cosmic_dbs_channel_to_type()
  } else if (sigclass == "SV") {
    vec_channel2type <- cosmic_sv_channel_to_type()
  } else if (sigclass == "RNA-SBS") {
    vec_channel2type <- cosmic_rna_sbs_channel_to_type()
  } else {
    stop("unexpected sigclass")
  }

  types <- vec_channel2type[match(channel, names(vec_channel2type))]
  return(types)
}

cosmic_sbs_channel_to_type <- function() {
  c(
    `A[C>A]A` = "C>A", `A[C>A]C` = "C>A", `A[C>A]G` = "C>A", `A[C>A]T` = "C>A",
    `A[C>G]A` = "C>G", `A[C>G]C` = "C>G", `A[C>G]G` = "C>G", `A[C>G]T` = "C>G",
    `A[C>T]A` = "C>T", `A[C>T]C` = "C>T", `A[C>T]G` = "C>T", `A[C>T]T` = "C>T",
    `A[T>A]A` = "T>A", `A[T>A]C` = "T>A", `A[T>A]G` = "T>A", `A[T>A]T` = "T>A",
    `A[T>C]A` = "T>C", `A[T>C]C` = "T>C", `A[T>C]G` = "T>C", `A[T>C]T` = "T>C",
    `A[T>G]A` = "T>G", `A[T>G]C` = "T>G", `A[T>G]G` = "T>G", `A[T>G]T` = "T>G",
    `C[C>A]A` = "C>A", `C[C>A]C` = "C>A", `C[C>A]G` = "C>A", `C[C>A]T` = "C>A",
    `C[C>G]A` = "C>G", `C[C>G]C` = "C>G", `C[C>G]G` = "C>G", `C[C>G]T` = "C>G",
    `C[C>T]A` = "C>T", `C[C>T]C` = "C>T", `C[C>T]G` = "C>T", `C[C>T]T` = "C>T",
    `C[T>A]A` = "T>A", `C[T>A]C` = "T>A", `C[T>A]G` = "T>A", `C[T>A]T` = "T>A",
    `C[T>C]A` = "T>C", `C[T>C]C` = "T>C", `C[T>C]G` = "T>C", `C[T>C]T` = "T>C",
    `C[T>G]A` = "T>G", `C[T>G]C` = "T>G", `C[T>G]G` = "T>G", `C[T>G]T` = "T>G",
    `G[C>A]A` = "C>A", `G[C>A]C` = "C>A", `G[C>A]G` = "C>A", `G[C>A]T` = "C>A",
    `G[C>G]A` = "C>G", `G[C>G]C` = "C>G", `G[C>G]G` = "C>G", `G[C>G]T` = "C>G",
    `G[C>T]A` = "C>T", `G[C>T]C` = "C>T", `G[C>T]G` = "C>T", `G[C>T]T` = "C>T",
    `G[T>A]A` = "T>A", `G[T>A]C` = "T>A", `G[T>A]G` = "T>A", `G[T>A]T` = "T>A",
    `G[T>C]A` = "T>C", `G[T>C]C` = "T>C", `G[T>C]G` = "T>C", `G[T>C]T` = "T>C",
    `G[T>G]A` = "T>G", `G[T>G]C` = "T>G", `G[T>G]G` = "T>G", `G[T>G]T` = "T>G",
    `T[C>A]A` = "C>A", `T[C>A]C` = "C>A", `T[C>A]G` = "C>A", `T[C>A]T` = "C>A",
    `T[C>G]A` = "C>G", `T[C>G]C` = "C>G", `T[C>G]G` = "C>G", `T[C>G]T` = "C>G",
    `T[C>T]A` = "C>T", `T[C>T]C` = "C>T", `T[C>T]G` = "C>T", `T[C>T]T` = "C>T",
    `T[T>A]A` = "T>A", `T[T>A]C` = "T>A", `T[T>A]G` = "T>A", `T[T>A]T` = "T>A",
    `T[T>C]A` = "T>C", `T[T>C]C` = "T>C", `T[T>C]G` = "T>C", `T[T>C]T` = "T>C",
    `T[T>G]A` = "T>G", `T[T>G]C` = "T>G", `T[T>G]G` = "T>G", `T[T>G]T` = "T>G"
  )
}


cosmic_cn_channel_to_type <- function() {
  c(
    `0:homdel:0-100kb` = "0", `0:homdel:100kb-1Mb` = "0", `0:homdel:>1Mb` = "0",
    `1:LOH:0-100kb` = "1", `1:LOH:100kb-1Mb` = "1", `1:LOH:1Mb-10Mb` = "1",
    `1:LOH:10Mb-40Mb` = "1", `1:LOH:>40Mb` = "1", `2:LOH:0-100kb` = "2",
    `2:LOH:100kb-1Mb` = "2", `2:LOH:1Mb-10Mb` = "2", `2:LOH:10Mb-40Mb` = "2",
    `2:LOH:>40Mb` = "2", `3-4:LOH:0-100kb` = "3-4", `3-4:LOH:100kb-1Mb` = "3-4",
    `3-4:LOH:1Mb-10Mb` = "3-4", `3-4:LOH:10Mb-40Mb` = "3-4", `3-4:LOH:>40Mb` = "3-4",
    `5-8:LOH:0-100kb` = "5-8", `5-8:LOH:100kb-1Mb` = "5-8", `5-8:LOH:1Mb-10Mb` = "5-8",
    `5-8:LOH:10Mb-40Mb` = "5-8", `5-8:LOH:>40Mb` = "5-8", `9+:LOH:0-100kb` = "9+",
    `9+:LOH:100kb-1Mb` = "9+", `9+:LOH:1Mb-10Mb` = "9+", `9+:LOH:10Mb-40Mb` = "9+",
    `9+:LOH:>40Mb` = "9+", `2:het:0-100kb` = "2", `2:het:100kb-1Mb` = "2",
    `2:het:1Mb-10Mb` = "2", `2:het:10Mb-40Mb` = "2", `2:het:>40Mb` = "2",
    `3-4:het:0-100kb` = "3-4", `3-4:het:100kb-1Mb` = "3-4", `3-4:het:1Mb-10Mb` = "3-4",
    `3-4:het:10Mb-40Mb` = "3-4", `3-4:het:>40Mb` = "3-4", `5-8:het:0-100kb` = "5-8",
    `5-8:het:100kb-1Mb` = "5-8", `5-8:het:1Mb-10Mb` = "5-8", `5-8:het:10Mb-40Mb` = "5-8",
    `5-8:het:>40Mb` = "5-8", `9+:het:0-100kb` = "9+", `9+:het:100kb-1Mb` = "9+",
    `9+:het:1Mb-10Mb` = "9+", `9+:het:10Mb-40Mb` = "9+", `9+:het:>40Mb` = "9+"
  )
}

cosmic_dbs_channel_to_type <- function() {
  c(
    `AC>CA` = "AC>NN", `AC>CG` = "AC>NN", `AC>CT` = "AC>NN", `AC>GA` = "AC>NN",
    `AC>GG` = "AC>NN", `AC>GT` = "AC>NN", `AC>TA` = "AC>NN", `AC>TG` = "AC>NN",
    `AC>TT` = "AC>NN", `AT>CA` = "AT>NN", `AT>CC` = "AT>NN", `AT>CG` = "AT>NN",
    `AT>GA` = "AT>NN", `AT>GC` = "AT>NN", `AT>TA` = "AT>NN", `CC>AA` = "CC>NN",
    `CC>AG` = "CC>NN", `CC>AT` = "CC>NN", `CC>GA` = "CC>NN", `CC>GG` = "CC>NN",
    `CC>GT` = "CC>NN", `CC>TA` = "CC>NN", `CC>TG` = "CC>NN", `CC>TT` = "CC>NN",
    `CG>AT` = "CG>NN", `CG>GC` = "CG>NN", `CG>GT` = "CG>NN", `CG>TA` = "CG>NN",
    `CG>TC` = "CG>NN", `CG>TT` = "CG>NN", `CT>AA` = "CT>NN", `CT>AC` = "CT>NN",
    `CT>AG` = "CT>NN", `CT>GA` = "CT>NN", `CT>GC` = "CT>NN", `CT>GG` = "CT>NN",
    `CT>TA` = "CT>NN", `CT>TC` = "CT>NN", `CT>TG` = "CT>NN", `GC>AA` = "GC>NN",
    `GC>AG` = "GC>NN", `GC>AT` = "GC>NN", `GC>CA` = "GC>NN", `GC>CG` = "GC>NN",
    `GC>TA` = "GC>NN", `TA>AT` = "TA>NN", `TA>CG` = "TA>NN", `TA>CT` = "TA>NN",
    `TA>GC` = "TA>NN", `TA>GG` = "TA>NN", `TA>GT` = "TA>NN", `TC>AA` = "TC>NN",
    `TC>AG` = "TC>NN", `TC>AT` = "TC>NN", `TC>CA` = "TC>NN", `TC>CG` = "TC>NN",
    `TC>CT` = "TC>NN", `TC>GA` = "TC>NN", `TC>GG` = "TC>NN", `TC>GT` = "TC>NN",
    `TG>AA` = "TG>NN", `TG>AC` = "TG>NN", `TG>AT` = "TG>NN", `TG>CA` = "TG>NN",
    `TG>CC` = "TG>NN", `TG>CT` = "TG>NN", `TG>GA` = "TG>NN", `TG>GC` = "TG>NN",
    `TG>GT` = "TG>NN", `TT>AA` = "TT>NN", `TT>AC` = "TT>NN", `TT>AG` = "TT>NN",
    `TT>CA` = "TT>NN", `TT>CC` = "TT>NN", `TT>CG` = "TT>NN", `TT>GA` = "TT>NN",
    `TT>GC` = "TT>NN", `TT>GG` = "TT>NN"
  )
}

cosmic_id_channel_to_type <- function() {
  c(
    `1:Del:C:0` = "1:Del:C", `1:Del:C:1` = "1:Del:C", `1:Del:C:2` = "1:Del:C",
    `1:Del:C:3` = "1:Del:C", `1:Del:C:4` = "1:Del:C", `1:Del:C:5` = "1:Del:C",
    `1:Del:T:0` = "1:Del:T", `1:Del:T:1` = "1:Del:T", `1:Del:T:2` = "1:Del:T",
    `1:Del:T:3` = "1:Del:T", `1:Del:T:4` = "1:Del:T", `1:Del:T:5` = "1:Del:T",
    `1:Ins:C:0` = "1:Ins:C", `1:Ins:C:1` = "1:Ins:C", `1:Ins:C:2` = "1:Ins:C",
    `1:Ins:C:3` = "1:Ins:C", `1:Ins:C:4` = "1:Ins:C", `1:Ins:C:5` = "1:Ins:C",
    `1:Ins:T:0` = "1:Ins:T", `1:Ins:T:1` = "1:Ins:T", `1:Ins:T:2` = "1:Ins:T",
    `1:Ins:T:3` = "1:Ins:T", `1:Ins:T:4` = "1:Ins:T", `1:Ins:T:5` = "1:Ins:T",
    `2:Del:R:0` = "2:Del:R", `2:Del:R:1` = "2:Del:R", `2:Del:R:2` = "2:Del:R",
    `2:Del:R:3` = "2:Del:R", `2:Del:R:4` = "2:Del:R", `2:Del:R:5` = "2:Del:R",
    `3:Del:R:0` = "3:Del:R", `3:Del:R:1` = "3:Del:R", `3:Del:R:2` = "3:Del:R",
    `3:Del:R:3` = "3:Del:R", `3:Del:R:4` = "3:Del:R", `3:Del:R:5` = "3:Del:R",
    `4:Del:R:0` = "4:Del:R", `4:Del:R:1` = "4:Del:R", `4:Del:R:2` = "4:Del:R",
    `4:Del:R:3` = "4:Del:R", `4:Del:R:4` = "4:Del:R", `4:Del:R:5` = "4:Del:R",
    `5:Del:R:0` = "5:Del:R", `5:Del:R:1` = "5:Del:R", `5:Del:R:2` = "5:Del:R",
    `5:Del:R:3` = "5:Del:R", `5:Del:R:4` = "5:Del:R", `5:Del:R:5` = "5:Del:R",
    `2:Ins:R:0` = "2:Ins:R", `2:Ins:R:1` = "2:Ins:R", `2:Ins:R:2` = "2:Ins:R",
    `2:Ins:R:3` = "2:Ins:R", `2:Ins:R:4` = "2:Ins:R", `2:Ins:R:5` = "2:Ins:R",
    `3:Ins:R:0` = "3:Ins:R", `3:Ins:R:1` = "3:Ins:R", `3:Ins:R:2` = "3:Ins:R",
    `3:Ins:R:3` = "3:Ins:R", `3:Ins:R:4` = "3:Ins:R", `3:Ins:R:5` = "3:Ins:R",
    `4:Ins:R:0` = "4:Ins:R", `4:Ins:R:1` = "4:Ins:R", `4:Ins:R:2` = "4:Ins:R",
    `4:Ins:R:3` = "4:Ins:R", `4:Ins:R:4` = "4:Ins:R", `4:Ins:R:5` = "4:Ins:R",
    `5:Ins:R:0` = "5:Ins:R", `5:Ins:R:1` = "5:Ins:R", `5:Ins:R:2` = "5:Ins:R",
    `5:Ins:R:3` = "5:Ins:R", `5:Ins:R:4` = "5:Ins:R", `5:Ins:R:5` = "5:Ins:R",
    `2:Del:M:1` = "2:Del:M", `3:Del:M:1` = "3:Del:M", `3:Del:M:2` = "3:Del:M",
    `4:Del:M:1` = "4:Del:M", `4:Del:M:2` = "4:Del:M", `4:Del:M:3` = "4:Del:M",
    `5:Del:M:1` = "5:Del:M", `5:Del:M:2` = "5:Del:M", `5:Del:M:3` = "5:Del:M",
    `5:Del:M:4` = "5:Del:M", `5:Del:M:5` = "5:Del:M"
  )
}

cosmic_sv_channel_to_type <- function(){
  c(
    "clustered_del_1-10Kb" = "clustered",
    "clustered_del_10-100Kb" = "clustered",
    "clustered_del_100Kb-1Mb" = "clustered",
    "clustered_del_1Mb-10Mb" = "clustered",
    "clustered_del_>10Mb" = "clustered",
    "clustered_tds_1-10Kb" = "clustered",
    "clustered_tds_10-100Kb" = "clustered",
    "clustered_tds_100Kb-1Mb" = "clustered",
    "clustered_tds_1Mb-10Mb" = "clustered",
    "clustered_tds_>10Mb" = "clustered",
    "clustered_inv_1-10Kb" = "clustered",
    "clustered_inv_10-100Kb" = "clustered",
    "clustered_inv_100Kb-1Mb" = "clustered",
    "clustered_inv_1Mb-10Mb" = "clustered",
    "clustered_inv_>10Mb" = "clustered",
    "clustered_trans" = "clustered",
    "non-clustered_del_1-10Kb" = "non-clustered",
    "non-clustered_del_10-100Kb" = "non-clustered",
    "non-clustered_del_100Kb-1Mb" = "non-clustered",
    "non-clustered_del_1Mb-10Mb" = "non-clustered",
    "non-clustered_del_>10Mb" = "non-clustered",
    "non-clustered_tds_1-10Kb" = "non-clustered",
    "non-clustered_tds_10-100Kb" = "non-clustered",
    "non-clustered_tds_100Kb-1Mb" = "non-clustered",
    "non-clustered_tds_1Mb-10Mb" = "non-clustered",
    "non-clustered_tds_>10Mb" = "non-clustered",
    "non-clustered_inv_1-10Kb" = "non-clustered",
    "non-clustered_inv_10-100Kb" = "non-clustered",
    "non-clustered_inv_100Kb-1Mb" = "non-clustered",
    "non-clustered_inv_1Mb-10Mb" = "non-clustered",
    "non-clustered_inv_>10Mb" = "non-clustered",
    "non-clustered_trans" = "non-clustered"
  )
}

cosmic_rna_sbs_channel_to_type <- function(){
  c(
    `A[A>C]A` = "A>C", `A[A>C]C` = "A>C", `A[A>C]G` = "A>C", `A[A>C]T` = "A>C",
    `A[A>G]A` = "A>G", `A[A>G]C` = "A>G", `A[A>G]G` = "A>G", `A[A>G]T` = "A>G",
    `A[A>T]A` = "A>T", `A[A>T]C` = "A>T", `A[A>T]G` = "A>T", `A[A>T]T` = "A>T",
    `A[C>A]A` = "C>A", `A[C>A]C` = "C>A", `A[C>A]G` = "C>A", `A[C>A]T` = "C>A",
    `A[C>G]A` = "C>G", `A[C>G]C` = "C>G", `A[C>G]G` = "C>G", `A[C>G]T` = "C>G",
    `A[C>T]A` = "C>T", `A[C>T]C` = "C>T", `A[C>T]G` = "C>T", `A[C>T]T` = "C>T",
    `A[G>A]A` = "G>A", `A[G>A]C` = "G>A", `A[G>A]G` = "G>A", `A[G>A]T` = "G>A",
    `A[G>C]A` = "G>C", `A[G>C]C` = "G>C", `A[G>C]G` = "G>C", `A[G>C]T` = "G>C",
    `A[G>T]A` = "G>T", `A[G>T]C` = "G>T", `A[G>T]G` = "G>T", `A[G>T]T` = "G>T",
    `A[T>A]A` = "T>A", `A[T>A]C` = "T>A", `A[T>A]G` = "T>A", `A[T>A]T` = "T>A",
    `A[T>C]A` = "T>C", `A[T>C]C` = "T>C", `A[T>C]G` = "T>C", `A[T>C]T` = "T>C",
    `A[T>G]A` = "T>G", `A[T>G]C` = "T>G", `A[T>G]G` = "T>G", `A[T>G]T` = "T>G",
    `C[A>C]A` = "A>C", `C[A>C]C` = "A>C", `C[A>C]G` = "A>C", `C[A>C]T` = "A>C",
    `C[A>G]A` = "A>G", `C[A>G]C` = "A>G", `C[A>G]G` = "A>G", `C[A>G]T` = "A>G",
    `C[A>T]A` = "A>T", `C[A>T]C` = "A>T", `C[A>T]G` = "A>T", `C[A>T]T` = "A>T",
    `C[C>A]A` = "C>A", `C[C>A]C` = "C>A", `C[C>A]G` = "C>A", `C[C>A]T` = "C>A",
    `C[C>G]A` = "C>G", `C[C>G]C` = "C>G", `C[C>G]G` = "C>G", `C[C>G]T` = "C>G",
    `C[C>T]A` = "C>T", `C[C>T]C` = "C>T", `C[C>T]G` = "C>T", `C[C>T]T` = "C>T",
    `C[G>A]A` = "G>A", `C[G>A]C` = "G>A", `C[G>A]G` = "G>A", `C[G>A]T` = "G>A",
    `C[G>C]A` = "G>C", `C[G>C]C` = "G>C", `C[G>C]G` = "G>C", `C[G>C]T` = "G>C",
    `C[G>T]A` = "G>T", `C[G>T]C` = "G>T", `C[G>T]G` = "G>T", `C[G>T]T` = "G>T",
    `C[T>A]A` = "T>A", `C[T>A]C` = "T>A", `C[T>A]G` = "T>A", `C[T>A]T` = "T>A",
    `C[T>C]A` = "T>C", `C[T>C]C` = "T>C", `C[T>C]G` = "T>C", `C[T>C]T` = "T>C",
    `C[T>G]A` = "T>G", `C[T>G]C` = "T>G", `C[T>G]G` = "T>G", `C[T>G]T` = "T>G",
    `G[A>C]A` = "A>C", `G[A>C]C` = "A>C", `G[A>C]G` = "A>C", `G[A>C]T` = "A>C",
    `G[A>G]A` = "A>G", `G[A>G]C` = "A>G", `G[A>G]G` = "A>G", `G[A>G]T` = "A>G",
    `G[A>T]A` = "A>T", `G[A>T]C` = "A>T", `G[A>T]G` = "A>T", `G[A>T]T` = "A>T",
    `G[C>A]A` = "C>A", `G[C>A]C` = "C>A", `G[C>A]G` = "C>A", `G[C>A]T` = "C>A",
    `G[C>G]A` = "C>G", `G[C>G]C` = "C>G", `G[C>G]G` = "C>G", `G[C>G]T` = "C>G",
    `G[C>T]A` = "C>T", `G[C>T]C` = "C>T", `G[C>T]G` = "C>T", `G[C>T]T` = "C>T",
    `G[G>A]A` = "G>A", `G[G>A]C` = "G>A", `G[G>A]G` = "G>A", `G[G>A]T` = "G>A",
    `G[G>C]A` = "G>C", `G[G>C]C` = "G>C", `G[G>C]G` = "G>C", `G[G>C]T` = "G>C",
    `G[G>T]A` = "G>T", `G[G>T]C` = "G>T", `G[G>T]G` = "G>T", `G[G>T]T` = "G>T",
    `G[T>A]A` = "T>A", `G[T>A]C` = "T>A", `G[T>A]G` = "T>A", `G[T>A]T` = "T>A",
    `G[T>C]A` = "T>C", `G[T>C]C` = "T>C", `G[T>C]G` = "T>C", `G[T>C]T` = "T>C",
    `G[T>G]A` = "T>G", `G[T>G]C` = "T>G", `G[T>G]G` = "T>G", `G[T>G]T` = "T>G",
    `T[A>C]A` = "A>C", `T[A>C]C` = "A>C", `T[A>C]G` = "A>C", `T[A>C]T` = "A>C",
    `T[A>G]A` = "A>G", `T[A>G]C` = "A>G", `T[A>G]G` = "A>G", `T[A>G]T` = "A>G",
    `T[A>T]A` = "A>T", `T[A>T]C` = "A>T", `T[A>T]G` = "A>T", `T[A>T]T` = "A>T",
    `T[C>A]A` = "C>A", `T[C>A]C` = "C>A", `T[C>A]G` = "C>A", `T[C>A]T` = "C>A",
    `T[C>G]A` = "C>G", `T[C>G]C` = "C>G", `T[C>G]G` = "C>G", `T[C>G]T` = "C>G",
    `T[C>T]A` = "C>T", `T[C>T]C` = "C>T", `T[C>T]G` = "C>T", `T[C>T]T` = "C>T",
    `T[G>A]A` = "G>A", `T[G>A]C` = "G>A", `T[G>A]G` = "G>A", `T[G>A]T` = "G>A",
    `T[G>C]A` = "G>C", `T[G>C]C` = "G>C", `T[G>C]G` = "G>C", `T[G>C]T` = "G>C",
    `T[G>T]A` = "G>T", `T[G>T]C` = "G>T", `T[G>T]G` = "G>T", `T[G>T]T` = "G>T",
    `T[T>A]A` = "T>A", `T[T>A]C` = "T>A", `T[T>A]G` = "T>A", `T[T>A]T` = "T>A",
    `T[T>C]A` = "T>C", `T[T>C]C` = "T>C", `T[T>C]G` = "T>C", `T[T>C]T` = "T>C",
    `T[T>G]A` = "T>G", `T[T>G]C` = "T>G", `T[T>G]G` = "T>G", `T[T>G]T` = "T>G"
  )
}
# Export Data -------------------------------------------------------------
#' Sigstash Collection -> Sigminer
#'
#' Converts a sigstash format signature collection to a sigminer-compatible signature database.
#' This lets you run sigminer signature analyses sigs using sigstash
#'
#' @param signatures a signature collection
#' @return a signature data.frame compatible with `sig` argument of sigmienr sig_fit()
#' @export
#'
#' @details
#' sigminer sig arguments allows a data.frame with rows
#' representing components and columns representing signatures
#'
#' @examples
#' # Load packages for signature fitting
#' library(sigminer)
#' library(quadprog)
#'
#' # Load available datasets
#' signature_collection <- sig_load("COSMIC_v3.3.1_SBS_GRCh38")
#'
#' # Sigminer
#' sigminer_collection <- sig_collection_to_sigminer(signature_collection)
#'
#' # Create Mock Sample Data
#' sample_data <- matrix(ncol = 3, runif(n = 96 * 3, min = 0, max = 20))
#'
#' # Run Sigminer
#' sig_fit(sample_data, sig = sigminer_collection)
#'
sig_collection_to_sigminer <- function(signatures) {
  sigshared::assert_signature_collection(signatures)
  assertions::assert_greater_than(length(signatures), minimum = 0)

  first_sig_channel_order <- signatures[[1]][["channel"]]

  ls <- lapply(seq_along(signatures), FUN = \(i){
    sig <- signatures[[i]]
    assertions::assert_identical(first_sig_channel_order, sig[["channel"]])
    df_fraction <- sig[, "fraction"]
    colnames(df_fraction) <- names(signatures)[i]
    return(df_fraction)
  })
  df_wide <- do.call("cbind", ls)

  # Capitalise kb to Kb in copynumber channel names
  first_sig_channel_order <- sub(x = first_sig_channel_order, "([0-9])kb", "\\1Kb")

  rownames(df_wide) <- first_sig_channel_order


  return(df_wide)
}
