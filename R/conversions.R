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
  types <- sig_convert_channel2type(channels, sigclass)

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

#' List Valid Signature Classes
#'
#' List sigverse supported signature feature-spaces.
#'
#' @return vector of sigverse supported signature catalogue feature-spaces.
#' @export
#'
#' @examples
#' sig_valid_sigclass()
sig_valid_sigclass <- function() {
  c("SBS96", "SBS1536", "ID83", "CN48", "CN80", "CN176", "DBS78", "SV32", "SV38", "RNA-SBS192")
}

#' Channel to Type
#'
#' Converts COSMIC-style channel names (e.g., "A\[C>A\]A") into higher-level mutation types (e.g., "C>A").
#'
#' @param channel A character vector of COSMIC-style feature channels (e.g., "A\[C>A\]A").
#' @param sigclass A string indicating the class of signatures these feature sets belong to. Run [sig_valid_sigclass()] to see all options.
#'
#' @return A character vector of the same length as `channel`, describing the higher-level mutation types.
#' @export
#'
#' @examples
#' doublet_channels <- c(
#'   "AC>CA", "AC>CG", "AC>CT", "AC>GA",
#'   "AC>GG", "AC>GT", "AC>TA", "AC>TG",
#'   "AC>TT", "AT>CA", "AT>CC", "AT>CG",
#'   "AT>GA", "AT>GC", "AT>TA", "CC>AA",
#'   "CC>AG", "CC>AT", "CC>GA", "CC>GG",
#'   "CC>GT", "CC>TA", "CC>TG", "CC>TT",
#'   "CG>AT", "CG>GC", "CG>GT", "CG>TA",
#'   "CG>TC", "CG>TT", "CT>AA", "CT>AC",
#'   "CT>AG", "CT>GA", "CT>GC", "CT>GG",
#'   "CT>TA", "CT>TC", "CT>TG", "GC>AA",
#'   "GC>AG", "GC>AT", "GC>CA", "GC>CG",
#'   "GC>TA", "TA>AT", "TA>CG", "TA>CT",
#'   "TA>GC", "TA>GG", "TA>GT", "TC>AA",
#'   "TC>AG", "TC>AT", "TC>CA", "TC>CG",
#'   "TC>CT", "TC>GA", "TC>GG", "TC>GT",
#'   "TG>AA", "TG>AC", "TG>AT", "TG>CA",
#'   "TG>CC", "TG>CT", "TG>GA", "TG>GC",
#'   "TG>GT", "TT>AA", "TT>AC", "TT>AG",
#'   "TT>CA", "TT>CC", "TT>CG", "TT>GA",
#'   "TT>GC", "TT>GG"
#' )
#'
#' # Get higher-level doublet channel types
#' sig_convert_channel2type(doublet_channels, sigclass = "DBS78")
#'
sig_convert_channel2type <- function(channel, sigclass = "SBS96") {
  # Assertions
  assertions::assert_string(sigclass)
  assertions::assert_subset(sigclass, sig_valid_sigclass())

  # Ensure channel is a character vector
  assertions::assert_character_vector(channel)

  # Retrieve the appropriate function for the given sigclass
  vec_channel2type <- sig_get_channel_to_type_maps(sigclass)

  # Map channels to their higher-level types using the retrieved function
  types <- vec_channel2type[match(channel, names(vec_channel2type))]
  unexpected_channels <- channel[is.na(types)]
  assertions::assert_no_missing(types, msg = "Failed to convert signature channels to 'types' due to unexpected channels for cosmic {sigclass} signatures: {unexpected_channels}. Run sig_get_valid_cosmic_channels('{sigclass}') to see a list of the expected channels")

  return(types)
}


#' Get Valid COSMIC Channels and Types
#'
#' Returns a character vector of all the expected channels or types for the specified signature class (sigclass).
#'
#' @param sigclass A string indicating the class of signatures. Run [sig_valid_sigclass()] to see all options.
#'
#' @return A character vector of valid COSMIC channels or types for the specified signature class.
#' @export
#'
#' @examples
#' sig_get_valid_cosmic_channels("SBS96")
#' sig_get_valid_cosmic_types("SBS96")
sig_get_valid_cosmic_channels <- function(sigclass) {
  # Assertions
  assertions::assert_string(sigclass)
  assertions::assert_subset(sigclass, sig_valid_sigclass())

  # Get Channels
  names(sig_get_channel_to_type_maps(sigclass))
}

#' @inherit sig_get_valid_cosmic_channels
#' @export
sig_get_valid_cosmic_types <- function(sigclass) {
  # Assertions
  assertions::assert_string(sigclass)
  assertions::assert_subset(sigclass, sig_valid_sigclass())

  unname(sig_get_channel_to_type_maps(sigclass))
}

sig_get_channel_to_type_maps <- function(sigclass) {
  vec_channel2type <- switch(sigclass,
    "SBS96" = cosmic_sbs96_channel_to_type(),
    "SBS1536" = cosmic_sbs1536_channel_to_type(),
    "ID83" = cosmic_id83_channel_to_type(),
    "CN48" = cosmic_cn48_channel_to_type(),
    "CN80" = cosmic_cn80_channel_to_type(),
    "CN176" = cosmic_cn176_channel_to_type(),
    "DBS78" = cosmic_dbs78_channel_to_type(),
    "SV32" = cosmic_sv32_channel_to_type(),
    "SV38" = cosmic_sv38_channel_to_type(),
    "RNA-SBS192" = cosmic_rna_sbs192_channel_to_type(),
    stop("unexpected sigclass")
  )
  return(vec_channel2type)
}


cosmic_sbs96_channel_to_type <- function() {
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


cosmic_cn48_channel_to_type <- function() {
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


cosmic_cn80_channel_to_type <- function() {
  c(`BP10MB[0]` = "BP10MB", `BP10MB[1]` = "BP10MB", `BP10MB[2]` = "BP10MB", `BP10MB[3]` = "BP10MB", `BP10MB[4]` = "BP10MB", `BP10MB[5]` = "BP10MB", `BP10MB[>5]` = "BP10MB", `BPArm[0]` = "BPArm", `BPArm[1]` = "BPArm", `BPArm[2]` = "BPArm", `BPArm[3]` = "BPArm", `BPArm[4]` = "BPArm", `BPArm[5]` = "BPArm", `BPArm[6]` = "BPArm", `BPArm[7]` = "BPArm", `BPArm[8]` = "BPArm", `BPArm[9]` = "BPArm", `BPArm[10]` = "BPArm", `BPArm[>10 & <=20]` = "BPArm", `BPArm[>20 & <=30]` = "BPArm", `BPArm[>30]` = "BPArm", `CN[0]` = "CN", `CN[1]` = "CN", `CN[2]` = "CN", `CN[3]` = "CN", `CN[4]` = "CN", `CN[>4 & <=8]` = "CN", `CN[>8]` = "CN", `CNCP[0]` = "CNCP", `CNCP[1]` = "CNCP", `CNCP[2]` = "CNCP", `CNCP[3]` = "CNCP", `CNCP[4]` = "CNCP", `CNCP[>4 & <=8]` = "CNCP", `CNCP[>8]` = "CNCP", `OsCN[0]` = "OsCN", `OsCN[1]` = "OsCN", `OsCN[2]` = "OsCN", `OsCN[3]` = "OsCN", `OsCN[4]` = "OsCN", `OsCN[>4 & <=10]` = "OsCN", `OsCN[>10]` = "OsCN", `SS[<=2]` = "SS", `SS[>2 & <=3]` = "SS", `SS[>3 & <=4]` = "SS", `SS[>4 & <=5]` = "SS", `SS[>5 & <=6]` = "SS", `SS[>6 & <=7]` = "SS", `SS[>7 & <=8]` = "SS", `SS[>8]` = "SS", `NC50[<=2]` = "NC50", `NC50[3]` = "NC50", `NC50[4]` = "NC50", `NC50[5]` = "NC50", `NC50[6]` = "NC50", `NC50[7]` = "NC50", `NC50[>7]` = "NC50", `BoChr[1]` = "BoChr", `BoChr[2]` = "BoChr", `BoChr[3]` = "BoChr", `BoChr[4]` = "BoChr", `BoChr[5]` = "BoChr", `BoChr[6]` = "BoChr", `BoChr[7]` = "BoChr", `BoChr[8]` = "BoChr", `BoChr[9]` = "BoChr", `BoChr[10]` = "BoChr", `BoChr[11]` = "BoChr", `BoChr[12]` = "BoChr", `BoChr[13]` = "BoChr", `BoChr[14]` = "BoChr", `BoChr[15]` = "BoChr", `BoChr[16]` = "BoChr", `BoChr[17]` = "BoChr", `BoChr[18]` = "BoChr", `BoChr[19]` = "BoChr", `BoChr[20]` = "BoChr", `BoChr[21]` = "BoChr", `BoChr[22]` = "BoChr", `BoChr[23]` = "BoChr")
}

cosmic_cn176_channel_to_type <- function() {
  c(`E:HH:0:AA` = "E:HH", `E:HH:0:BB` = "E:HH", `E:HH:1:AA` = "E:HH", `E:HH:1:BB` = "E:HH", `E:HH:2:AA` = "E:HH", `E:HH:2:BB` = "E:HH", `E:HH:3:AA` = "E:HH", `E:HH:3:BB` = "E:HH", `E:HH:4:AA` = "E:HH", `E:HH:4:BB` = "E:HH", `E:HH:5-8:AA` = "E:HH", `E:HH:5-8:BB` = "E:HH", `E:HH:9+:AA` = "E:HH", `E:HH:9+:BB` = "E:HH", `E:HL:1:AA` = "E:HL", `E:HL:2:AA` = "E:HL", `E:HL:3:AA` = "E:HL", `E:HL:3:BB` = "E:HL", `E:HL:4:AA` = "E:HL", `E:HL:4:BB` = "E:HL", `E:HL:5-8:AA` = "E:HL", `E:HL:5-8:BB` = "E:HL", `E:HL:9+:AA` = "E:HL", `E:HL:9+:BB` = "E:HL", `E:LH:1:AA` = "E:LH", `E:LH:2:AA` = "E:LH", `E:LH:3:AA` = "E:LH", `E:LH:3:BB` = "E:LH", `E:LH:4:AA` = "E:LH", `E:LH:4:BB` = "E:LH", `E:LH:5-8:AA` = "E:LH", `E:LH:5-8:BB` = "E:LH", `E:LH:9+:AA` = "E:LH", `E:LH:9+:BB` = "E:LH", `E:LL:1:AA` = "E:LL", `E:LL:2:AA` = "E:LL", `E:LL:3:AA` = "E:LL", `E:LL:3:BB` = "E:LL", `E:LL:4:AA` = "E:LL", `E:LL:4:BB` = "E:LL", `E:LL:5-8:AA` = "E:LL", `E:LL:5-8:BB` = "E:LL", `E:LL:9+:AA` = "E:LL", `E:LL:9+:BB` = "E:LL", `L:HH:0:AA` = "L:HH", `L:HH:0:BB` = "L:HH", `L:HH:1:AA` = "L:HH", `L:HH:1:BB` = "L:HH", `L:HH:2:AA` = "L:HH", `L:HH:2:BB` = "L:HH", `L:HH:3:AA` = "L:HH", `L:HH:3:BB` = "L:HH", `L:HH:4:AA` = "L:HH", `L:HH:4:BB` = "L:HH", `L:HH:5-8:AA` = "L:HH", `L:HH:5-8:BB` = "L:HH", `L:HH:9+:AA` = "L:HH", `L:HH:9+:BB` = "L:HH", `L:HL:1:AA` = "L:HL", `L:HL:2:AA` = "L:HL", `L:HL:3:AA` = "L:HL", `L:HL:3:BB` = "L:HL", `L:HL:4:AA` = "L:HL", `L:HL:4:BB` = "L:HL", `L:HL:5-8:AA` = "L:HL", `L:HL:5-8:BB` = "L:HL", `L:HL:9+:AA` = "L:HL", `L:HL:9+:BB` = "L:HL", `L:LH:1:AA` = "L:LH", `L:LH:2:AA` = "L:LH", `L:LH:3:AA` = "L:LH", `L:LH:3:BB` = "L:LH", `L:LH:4:AA` = "L:LH", `L:LH:4:BB` = "L:LH", `L:LH:5-8:AA` = "L:LH", `L:LH:5-8:BB` = "L:LH", `L:LH:9+:AA` = "L:LH", `L:LH:9+:BB` = "L:LH", `L:LL:1:AA` = "L:LL", `L:LL:2:AA` = "L:LL", `L:LL:3:AA` = "L:LL", `L:LL:3:BB` = "L:LL", `L:LL:4:AA` = "L:LL", `L:LL:4:BB` = "L:LL", `L:LL:5-8:AA` = "L:LL", `L:LL:5-8:BB` = "L:LL", `L:LL:9+:AA` = "L:LL", `L:LL:9+:BB` = "L:LL", `M:HH:0:AA` = "M:HH", `M:HH:0:BB` = "M:HH", `M:HH:1:AA` = "M:HH", `M:HH:1:BB` = "M:HH", `M:HH:2:AA` = "M:HH", `M:HH:2:BB` = "M:HH", `M:HH:3:AA` = "M:HH", `M:HH:3:BB` = "M:HH", `M:HH:4:AA` = "M:HH", `M:HH:4:BB` = "M:HH", `M:HH:5-8:AA` = "M:HH", `M:HH:5-8:BB` = "M:HH", `M:HH:9+:AA` = "M:HH", `M:HH:9+:BB` = "M:HH", `M:HL:1:AA` = "M:HL", `M:HL:2:AA` = "M:HL", `M:HL:3:AA` = "M:HL", `M:HL:3:BB` = "M:HL", `M:HL:4:AA` = "M:HL", `M:HL:4:BB` = "M:HL", `M:HL:5-8:AA` = "M:HL", `M:HL:5-8:BB` = "M:HL", `M:HL:9+:AA` = "M:HL", `M:HL:9+:BB` = "M:HL", `M:LH:1:AA` = "M:LH", `M:LH:2:AA` = "M:LH", `M:LH:3:AA` = "M:LH", `M:LH:3:BB` = "M:LH", `M:LH:4:AA` = "M:LH", `M:LH:4:BB` = "M:LH", `M:LH:5-8:AA` = "M:LH", `M:LH:5-8:BB` = "M:LH", `M:LH:9+:AA` = "M:LH", `M:LH:9+:BB` = "M:LH", `M:LL:1:AA` = "M:LL", `M:LL:2:AA` = "M:LL", `M:LL:3:AA` = "M:LL", `M:LL:3:BB` = "M:LL", `M:LL:4:AA` = "M:LL", `M:LL:4:BB` = "M:LL", `M:LL:5-8:AA` = "M:LL", `M:LL:5-8:BB` = "M:LL", `M:LL:9+:AA` = "M:LL", `M:LL:9+:BB` = "M:LL", `S:HH:0:AA` = "S:HH", `S:HH:0:BB` = "S:HH", `S:HH:1:AA` = "S:HH", `S:HH:1:BB` = "S:HH", `S:HH:2:AA` = "S:HH", `S:HH:2:BB` = "S:HH", `S:HH:3:AA` = "S:HH", `S:HH:3:BB` = "S:HH", `S:HH:4:AA` = "S:HH", `S:HH:4:BB` = "S:HH", `S:HH:5-8:AA` = "S:HH", `S:HH:5-8:BB` = "S:HH", `S:HH:9+:AA` = "S:HH", `S:HH:9+:BB` = "S:HH", `S:HL:1:AA` = "S:HL", `S:HL:2:AA` = "S:HL", `S:HL:3:AA` = "S:HL", `S:HL:3:BB` = "S:HL", `S:HL:4:AA` = "S:HL", `S:HL:4:BB` = "S:HL", `S:HL:5-8:AA` = "S:HL", `S:HL:5-8:BB` = "S:HL", `S:HL:9+:AA` = "S:HL", `S:HL:9+:BB` = "S:HL", `S:LH:1:AA` = "S:LH", `S:LH:2:AA` = "S:LH", `S:LH:3:AA` = "S:LH", `S:LH:3:BB` = "S:LH", `S:LH:4:AA` = "S:LH", `S:LH:4:BB` = "S:LH", `S:LH:5-8:AA` = "S:LH", `S:LH:5-8:BB` = "S:LH", `S:LH:9+:AA` = "S:LH", `S:LH:9+:BB` = "S:LH", `S:LL:1:AA` = "S:LL", `S:LL:2:AA` = "S:LL", `S:LL:3:AA` = "S:LL", `S:LL:3:BB` = "S:LL", `S:LL:4:AA` = "S:LL", `S:LL:4:BB` = "S:LL", `S:LL:5-8:AA` = "S:LL", `S:LL:5-8:BB` = "S:LL", `S:LL:9+:AA` = "S:LL", `S:LL:9+:BB` = "S:LL")
}

cosmic_sbs1536_channel_to_type <- function() {
  c(`AA[C>A]AA` = "C>A", `AA[C>A]AC` = "C>A", `AA[C>A]AG` = "C>A", `AA[C>A]AT` = "C>A", `AA[C>A]CA` = "C>A", `AA[C>A]CC` = "C>A", `AA[C>A]CG` = "C>A", `AA[C>A]CT` = "C>A", `AA[C>A]GA` = "C>A", `AA[C>A]GC` = "C>A", `AA[C>A]GG` = "C>A", `AA[C>A]GT` = "C>A", `AA[C>A]TA` = "C>A", `AA[C>A]TC` = "C>A", `AA[C>A]TG` = "C>A", `AA[C>A]TT` = "C>A", `AA[C>G]AA` = "C>G", `AA[C>G]AC` = "C>G", `AA[C>G]AG` = "C>G", `AA[C>G]AT` = "C>G", `AA[C>G]CA` = "C>G", `AA[C>G]CC` = "C>G", `AA[C>G]CG` = "C>G", `AA[C>G]CT` = "C>G", `AA[C>G]GA` = "C>G", `AA[C>G]GC` = "C>G", `AA[C>G]GG` = "C>G", `AA[C>G]GT` = "C>G", `AA[C>G]TA` = "C>G", `AA[C>G]TC` = "C>G", `AA[C>G]TG` = "C>G", `AA[C>G]TT` = "C>G", `AA[C>T]AA` = "C>T", `AA[C>T]AC` = "C>T", `AA[C>T]AG` = "C>T", `AA[C>T]AT` = "C>T", `AA[C>T]CA` = "C>T", `AA[C>T]CC` = "C>T", `AA[C>T]CG` = "C>T", `AA[C>T]CT` = "C>T", `AA[C>T]GA` = "C>T", `AA[C>T]GC` = "C>T", `AA[C>T]GG` = "C>T", `AA[C>T]GT` = "C>T", `AA[C>T]TA` = "C>T", `AA[C>T]TC` = "C>T", `AA[C>T]TG` = "C>T", `AA[C>T]TT` = "C>T", `AA[T>A]AA` = "T>A", `AA[T>A]AC` = "T>A", `AA[T>A]AG` = "T>A", `AA[T>A]AT` = "T>A", `AA[T>A]CA` = "T>A", `AA[T>A]CC` = "T>A", `AA[T>A]CG` = "T>A", `AA[T>A]CT` = "T>A", `AA[T>A]GA` = "T>A", `AA[T>A]GC` = "T>A", `AA[T>A]GG` = "T>A", `AA[T>A]GT` = "T>A", `AA[T>A]TA` = "T>A", `AA[T>A]TC` = "T>A", `AA[T>A]TG` = "T>A", `AA[T>A]TT` = "T>A", `AA[T>C]AA` = "T>C", `AA[T>C]AC` = "T>C", `AA[T>C]AG` = "T>C", `AA[T>C]AT` = "T>C", `AA[T>C]CA` = "T>C", `AA[T>C]CC` = "T>C", `AA[T>C]CG` = "T>C", `AA[T>C]CT` = "T>C", `AA[T>C]GA` = "T>C", `AA[T>C]GC` = "T>C", `AA[T>C]GG` = "T>C", `AA[T>C]GT` = "T>C", `AA[T>C]TA` = "T>C", `AA[T>C]TC` = "T>C", `AA[T>C]TG` = "T>C", `AA[T>C]TT` = "T>C", `AA[T>G]AA` = "T>G", `AA[T>G]AC` = "T>G", `AA[T>G]AG` = "T>G", `AA[T>G]AT` = "T>G", `AA[T>G]CA` = "T>G", `AA[T>G]CC` = "T>G", `AA[T>G]CG` = "T>G", `AA[T>G]CT` = "T>G", `AA[T>G]GA` = "T>G", `AA[T>G]GC` = "T>G", `AA[T>G]GG` = "T>G", `AA[T>G]GT` = "T>G", `AA[T>G]TA` = "T>G", `AA[T>G]TC` = "T>G", `AA[T>G]TG` = "T>G", `AA[T>G]TT` = "T>G", `AC[C>A]AA` = "C>A", `AC[C>A]AC` = "C>A", `AC[C>A]AG` = "C>A", `AC[C>A]AT` = "C>A", `AC[C>A]CA` = "C>A", `AC[C>A]CC` = "C>A", `AC[C>A]CG` = "C>A", `AC[C>A]CT` = "C>A", `AC[C>A]GA` = "C>A", `AC[C>A]GC` = "C>A", `AC[C>A]GG` = "C>A", `AC[C>A]GT` = "C>A", `AC[C>A]TA` = "C>A", `AC[C>A]TC` = "C>A", `AC[C>A]TG` = "C>A", `AC[C>A]TT` = "C>A", `AC[C>G]AA` = "C>G", `AC[C>G]AC` = "C>G", `AC[C>G]AG` = "C>G", `AC[C>G]AT` = "C>G", `AC[C>G]CA` = "C>G", `AC[C>G]CC` = "C>G", `AC[C>G]CG` = "C>G", `AC[C>G]CT` = "C>G", `AC[C>G]GA` = "C>G", `AC[C>G]GC` = "C>G", `AC[C>G]GG` = "C>G", `AC[C>G]GT` = "C>G", `AC[C>G]TA` = "C>G", `AC[C>G]TC` = "C>G", `AC[C>G]TG` = "C>G", `AC[C>G]TT` = "C>G", `AC[C>T]AA` = "C>T", `AC[C>T]AC` = "C>T", `AC[C>T]AG` = "C>T", `AC[C>T]AT` = "C>T", `AC[C>T]CA` = "C>T", `AC[C>T]CC` = "C>T", `AC[C>T]CG` = "C>T", `AC[C>T]CT` = "C>T", `AC[C>T]GA` = "C>T", `AC[C>T]GC` = "C>T", `AC[C>T]GG` = "C>T", `AC[C>T]GT` = "C>T", `AC[C>T]TA` = "C>T", `AC[C>T]TC` = "C>T", `AC[C>T]TG` = "C>T", `AC[C>T]TT` = "C>T", `AC[T>A]AA` = "T>A", `AC[T>A]AC` = "T>A", `AC[T>A]AG` = "T>A", `AC[T>A]AT` = "T>A", `AC[T>A]CA` = "T>A", `AC[T>A]CC` = "T>A", `AC[T>A]CG` = "T>A", `AC[T>A]CT` = "T>A", `AC[T>A]GA` = "T>A", `AC[T>A]GC` = "T>A", `AC[T>A]GG` = "T>A", `AC[T>A]GT` = "T>A", `AC[T>A]TA` = "T>A", `AC[T>A]TC` = "T>A", `AC[T>A]TG` = "T>A", `AC[T>A]TT` = "T>A", `AC[T>C]AA` = "T>C", `AC[T>C]AC` = "T>C", `AC[T>C]AG` = "T>C", `AC[T>C]AT` = "T>C", `AC[T>C]CA` = "T>C", `AC[T>C]CC` = "T>C", `AC[T>C]CG` = "T>C", `AC[T>C]CT` = "T>C", `AC[T>C]GA` = "T>C", `AC[T>C]GC` = "T>C", `AC[T>C]GG` = "T>C", `AC[T>C]GT` = "T>C", `AC[T>C]TA` = "T>C", `AC[T>C]TC` = "T>C", `AC[T>C]TG` = "T>C", `AC[T>C]TT` = "T>C", `AC[T>G]AA` = "T>G", `AC[T>G]AC` = "T>G", `AC[T>G]AG` = "T>G", `AC[T>G]AT` = "T>G", `AC[T>G]CA` = "T>G", `AC[T>G]CC` = "T>G", `AC[T>G]CG` = "T>G", `AC[T>G]CT` = "T>G", `AC[T>G]GA` = "T>G", `AC[T>G]GC` = "T>G", `AC[T>G]GG` = "T>G", `AC[T>G]GT` = "T>G", `AC[T>G]TA` = "T>G", `AC[T>G]TC` = "T>G", `AC[T>G]TG` = "T>G", `AC[T>G]TT` = "T>G", `AG[C>A]AA` = "C>A", `AG[C>A]AC` = "C>A", `AG[C>A]AG` = "C>A", `AG[C>A]AT` = "C>A", `AG[C>A]CA` = "C>A", `AG[C>A]CC` = "C>A", `AG[C>A]CG` = "C>A", `AG[C>A]CT` = "C>A", `AG[C>A]GA` = "C>A", `AG[C>A]GC` = "C>A", `AG[C>A]GG` = "C>A", `AG[C>A]GT` = "C>A", `AG[C>A]TA` = "C>A", `AG[C>A]TC` = "C>A", `AG[C>A]TG` = "C>A", `AG[C>A]TT` = "C>A", `AG[C>G]AA` = "C>G", `AG[C>G]AC` = "C>G", `AG[C>G]AG` = "C>G", `AG[C>G]AT` = "C>G", `AG[C>G]CA` = "C>G", `AG[C>G]CC` = "C>G", `AG[C>G]CG` = "C>G", `AG[C>G]CT` = "C>G", `AG[C>G]GA` = "C>G", `AG[C>G]GC` = "C>G", `AG[C>G]GG` = "C>G", `AG[C>G]GT` = "C>G", `AG[C>G]TA` = "C>G", `AG[C>G]TC` = "C>G", `AG[C>G]TG` = "C>G", `AG[C>G]TT` = "C>G", `AG[C>T]AA` = "C>T", `AG[C>T]AC` = "C>T", `AG[C>T]AG` = "C>T", `AG[C>T]AT` = "C>T", `AG[C>T]CA` = "C>T", `AG[C>T]CC` = "C>T", `AG[C>T]CG` = "C>T", `AG[C>T]CT` = "C>T", `AG[C>T]GA` = "C>T", `AG[C>T]GC` = "C>T", `AG[C>T]GG` = "C>T", `AG[C>T]GT` = "C>T", `AG[C>T]TA` = "C>T", `AG[C>T]TC` = "C>T", `AG[C>T]TG` = "C>T", `AG[C>T]TT` = "C>T", `AG[T>A]AA` = "T>A", `AG[T>A]AC` = "T>A", `AG[T>A]AG` = "T>A", `AG[T>A]AT` = "T>A", `AG[T>A]CA` = "T>A", `AG[T>A]CC` = "T>A", `AG[T>A]CG` = "T>A", `AG[T>A]CT` = "T>A", `AG[T>A]GA` = "T>A", `AG[T>A]GC` = "T>A", `AG[T>A]GG` = "T>A", `AG[T>A]GT` = "T>A", `AG[T>A]TA` = "T>A", `AG[T>A]TC` = "T>A", `AG[T>A]TG` = "T>A", `AG[T>A]TT` = "T>A", `AG[T>C]AA` = "T>C", `AG[T>C]AC` = "T>C", `AG[T>C]AG` = "T>C", `AG[T>C]AT` = "T>C", `AG[T>C]CA` = "T>C", `AG[T>C]CC` = "T>C", `AG[T>C]CG` = "T>C", `AG[T>C]CT` = "T>C", `AG[T>C]GA` = "T>C", `AG[T>C]GC` = "T>C", `AG[T>C]GG` = "T>C", `AG[T>C]GT` = "T>C", `AG[T>C]TA` = "T>C", `AG[T>C]TC` = "T>C", `AG[T>C]TG` = "T>C", `AG[T>C]TT` = "T>C", `AG[T>G]AA` = "T>G", `AG[T>G]AC` = "T>G", `AG[T>G]AG` = "T>G", `AG[T>G]AT` = "T>G", `AG[T>G]CA` = "T>G", `AG[T>G]CC` = "T>G", `AG[T>G]CG` = "T>G", `AG[T>G]CT` = "T>G", `AG[T>G]GA` = "T>G", `AG[T>G]GC` = "T>G", `AG[T>G]GG` = "T>G", `AG[T>G]GT` = "T>G", `AG[T>G]TA` = "T>G", `AG[T>G]TC` = "T>G", `AG[T>G]TG` = "T>G", `AG[T>G]TT` = "T>G", `AT[C>A]AA` = "C>A", `AT[C>A]AC` = "C>A", `AT[C>A]AG` = "C>A", `AT[C>A]AT` = "C>A", `AT[C>A]CA` = "C>A", `AT[C>A]CC` = "C>A", `AT[C>A]CG` = "C>A", `AT[C>A]CT` = "C>A", `AT[C>A]GA` = "C>A", `AT[C>A]GC` = "C>A", `AT[C>A]GG` = "C>A", `AT[C>A]GT` = "C>A", `AT[C>A]TA` = "C>A", `AT[C>A]TC` = "C>A", `AT[C>A]TG` = "C>A", `AT[C>A]TT` = "C>A", `AT[C>G]AA` = "C>G", `AT[C>G]AC` = "C>G", `AT[C>G]AG` = "C>G", `AT[C>G]AT` = "C>G", `AT[C>G]CA` = "C>G", `AT[C>G]CC` = "C>G", `AT[C>G]CG` = "C>G", `AT[C>G]CT` = "C>G", `AT[C>G]GA` = "C>G", `AT[C>G]GC` = "C>G", `AT[C>G]GG` = "C>G", `AT[C>G]GT` = "C>G", `AT[C>G]TA` = "C>G", `AT[C>G]TC` = "C>G", `AT[C>G]TG` = "C>G", `AT[C>G]TT` = "C>G", `AT[C>T]AA` = "C>T", `AT[C>T]AC` = "C>T", `AT[C>T]AG` = "C>T", `AT[C>T]AT` = "C>T", `AT[C>T]CA` = "C>T", `AT[C>T]CC` = "C>T", `AT[C>T]CG` = "C>T", `AT[C>T]CT` = "C>T", `AT[C>T]GA` = "C>T", `AT[C>T]GC` = "C>T", `AT[C>T]GG` = "C>T", `AT[C>T]GT` = "C>T", `AT[C>T]TA` = "C>T", `AT[C>T]TC` = "C>T", `AT[C>T]TG` = "C>T", `AT[C>T]TT` = "C>T", `AT[T>A]AA` = "T>A", `AT[T>A]AC` = "T>A", `AT[T>A]AG` = "T>A", `AT[T>A]AT` = "T>A", `AT[T>A]CA` = "T>A", `AT[T>A]CC` = "T>A", `AT[T>A]CG` = "T>A", `AT[T>A]CT` = "T>A", `AT[T>A]GA` = "T>A", `AT[T>A]GC` = "T>A", `AT[T>A]GG` = "T>A", `AT[T>A]GT` = "T>A", `AT[T>A]TA` = "T>A", `AT[T>A]TC` = "T>A", `AT[T>A]TG` = "T>A", `AT[T>A]TT` = "T>A", `AT[T>C]AA` = "T>C", `AT[T>C]AC` = "T>C", `AT[T>C]AG` = "T>C", `AT[T>C]AT` = "T>C", `AT[T>C]CA` = "T>C", `AT[T>C]CC` = "T>C", `AT[T>C]CG` = "T>C", `AT[T>C]CT` = "T>C", `AT[T>C]GA` = "T>C", `AT[T>C]GC` = "T>C", `AT[T>C]GG` = "T>C", `AT[T>C]GT` = "T>C", `AT[T>C]TA` = "T>C", `AT[T>C]TC` = "T>C", `AT[T>C]TG` = "T>C", `AT[T>C]TT` = "T>C", `AT[T>G]AA` = "T>G", `AT[T>G]AC` = "T>G", `AT[T>G]AG` = "T>G", `AT[T>G]AT` = "T>G", `AT[T>G]CA` = "T>G", `AT[T>G]CC` = "T>G", `AT[T>G]CG` = "T>G", `AT[T>G]CT` = "T>G", `AT[T>G]GA` = "T>G", `AT[T>G]GC` = "T>G", `AT[T>G]GG` = "T>G", `AT[T>G]GT` = "T>G", `AT[T>G]TA` = "T>G", `AT[T>G]TC` = "T>G", `AT[T>G]TG` = "T>G", `AT[T>G]TT` = "T>G", `CA[C>A]AA` = "C>A", `CA[C>A]AC` = "C>A", `CA[C>A]AG` = "C>A", `CA[C>A]AT` = "C>A", `CA[C>A]CA` = "C>A", `CA[C>A]CC` = "C>A", `CA[C>A]CG` = "C>A", `CA[C>A]CT` = "C>A", `CA[C>A]GA` = "C>A", `CA[C>A]GC` = "C>A", `CA[C>A]GG` = "C>A", `CA[C>A]GT` = "C>A", `CA[C>A]TA` = "C>A", `CA[C>A]TC` = "C>A", `CA[C>A]TG` = "C>A", `CA[C>A]TT` = "C>A", `CA[C>G]AA` = "C>G", `CA[C>G]AC` = "C>G", `CA[C>G]AG` = "C>G", `CA[C>G]AT` = "C>G", `CA[C>G]CA` = "C>G", `CA[C>G]CC` = "C>G", `CA[C>G]CG` = "C>G", `CA[C>G]CT` = "C>G", `CA[C>G]GA` = "C>G", `CA[C>G]GC` = "C>G", `CA[C>G]GG` = "C>G", `CA[C>G]GT` = "C>G", `CA[C>G]TA` = "C>G", `CA[C>G]TC` = "C>G", `CA[C>G]TG` = "C>G", `CA[C>G]TT` = "C>G", `CA[C>T]AA` = "C>T", `CA[C>T]AC` = "C>T", `CA[C>T]AG` = "C>T", `CA[C>T]AT` = "C>T", `CA[C>T]CA` = "C>T", `CA[C>T]CC` = "C>T", `CA[C>T]CG` = "C>T", `CA[C>T]CT` = "C>T", `CA[C>T]GA` = "C>T", `CA[C>T]GC` = "C>T", `CA[C>T]GG` = "C>T", `CA[C>T]GT` = "C>T", `CA[C>T]TA` = "C>T", `CA[C>T]TC` = "C>T", `CA[C>T]TG` = "C>T", `CA[C>T]TT` = "C>T", `CA[T>A]AA` = "T>A", `CA[T>A]AC` = "T>A", `CA[T>A]AG` = "T>A", `CA[T>A]AT` = "T>A", `CA[T>A]CA` = "T>A", `CA[T>A]CC` = "T>A", `CA[T>A]CG` = "T>A", `CA[T>A]CT` = "T>A", `CA[T>A]GA` = "T>A", `CA[T>A]GC` = "T>A", `CA[T>A]GG` = "T>A", `CA[T>A]GT` = "T>A", `CA[T>A]TA` = "T>A", `CA[T>A]TC` = "T>A", `CA[T>A]TG` = "T>A", `CA[T>A]TT` = "T>A", `CA[T>C]AA` = "T>C", `CA[T>C]AC` = "T>C", `CA[T>C]AG` = "T>C", `CA[T>C]AT` = "T>C", `CA[T>C]CA` = "T>C", `CA[T>C]CC` = "T>C", `CA[T>C]CG` = "T>C", `CA[T>C]CT` = "T>C", `CA[T>C]GA` = "T>C", `CA[T>C]GC` = "T>C", `CA[T>C]GG` = "T>C", `CA[T>C]GT` = "T>C", `CA[T>C]TA` = "T>C", `CA[T>C]TC` = "T>C", `CA[T>C]TG` = "T>C", `CA[T>C]TT` = "T>C", `CA[T>G]AA` = "T>G", `CA[T>G]AC` = "T>G", `CA[T>G]AG` = "T>G", `CA[T>G]AT` = "T>G", `CA[T>G]CA` = "T>G", `CA[T>G]CC` = "T>G", `CA[T>G]CG` = "T>G", `CA[T>G]CT` = "T>G", `CA[T>G]GA` = "T>G", `CA[T>G]GC` = "T>G", `CA[T>G]GG` = "T>G", `CA[T>G]GT` = "T>G", `CA[T>G]TA` = "T>G", `CA[T>G]TC` = "T>G", `CA[T>G]TG` = "T>G", `CA[T>G]TT` = "T>G", `CC[C>A]AA` = "C>A", `CC[C>A]AC` = "C>A", `CC[C>A]AG` = "C>A", `CC[C>A]AT` = "C>A", `CC[C>A]CA` = "C>A", `CC[C>A]CC` = "C>A", `CC[C>A]CG` = "C>A", `CC[C>A]CT` = "C>A", `CC[C>A]GA` = "C>A", `CC[C>A]GC` = "C>A", `CC[C>A]GG` = "C>A", `CC[C>A]GT` = "C>A", `CC[C>A]TA` = "C>A", `CC[C>A]TC` = "C>A", `CC[C>A]TG` = "C>A", `CC[C>A]TT` = "C>A", `CC[C>G]AA` = "C>G", `CC[C>G]AC` = "C>G", `CC[C>G]AG` = "C>G", `CC[C>G]AT` = "C>G", `CC[C>G]CA` = "C>G", `CC[C>G]CC` = "C>G", `CC[C>G]CG` = "C>G", `CC[C>G]CT` = "C>G", `CC[C>G]GA` = "C>G", `CC[C>G]GC` = "C>G", `CC[C>G]GG` = "C>G", `CC[C>G]GT` = "C>G", `CC[C>G]TA` = "C>G", `CC[C>G]TC` = "C>G", `CC[C>G]TG` = "C>G", `CC[C>G]TT` = "C>G", `CC[C>T]AA` = "C>T", `CC[C>T]AC` = "C>T", `CC[C>T]AG` = "C>T", `CC[C>T]AT` = "C>T", `CC[C>T]CA` = "C>T", `CC[C>T]CC` = "C>T", `CC[C>T]CG` = "C>T", `CC[C>T]CT` = "C>T", `CC[C>T]GA` = "C>T", `CC[C>T]GC` = "C>T", `CC[C>T]GG` = "C>T", `CC[C>T]GT` = "C>T", `CC[C>T]TA` = "C>T", `CC[C>T]TC` = "C>T", `CC[C>T]TG` = "C>T", `CC[C>T]TT` = "C>T", `CC[T>A]AA` = "T>A", `CC[T>A]AC` = "T>A", `CC[T>A]AG` = "T>A", `CC[T>A]AT` = "T>A", `CC[T>A]CA` = "T>A", `CC[T>A]CC` = "T>A", `CC[T>A]CG` = "T>A", `CC[T>A]CT` = "T>A", `CC[T>A]GA` = "T>A", `CC[T>A]GC` = "T>A", `CC[T>A]GG` = "T>A", `CC[T>A]GT` = "T>A", `CC[T>A]TA` = "T>A", `CC[T>A]TC` = "T>A", `CC[T>A]TG` = "T>A", `CC[T>A]TT` = "T>A", `CC[T>C]AA` = "T>C", `CC[T>C]AC` = "T>C", `CC[T>C]AG` = "T>C", `CC[T>C]AT` = "T>C", `CC[T>C]CA` = "T>C", `CC[T>C]CC` = "T>C", `CC[T>C]CG` = "T>C", `CC[T>C]CT` = "T>C", `CC[T>C]GA` = "T>C", `CC[T>C]GC` = "T>C", `CC[T>C]GG` = "T>C", `CC[T>C]GT` = "T>C", `CC[T>C]TA` = "T>C", `CC[T>C]TC` = "T>C", `CC[T>C]TG` = "T>C", `CC[T>C]TT` = "T>C", `CC[T>G]AA` = "T>G", `CC[T>G]AC` = "T>G", `CC[T>G]AG` = "T>G", `CC[T>G]AT` = "T>G", `CC[T>G]CA` = "T>G", `CC[T>G]CC` = "T>G", `CC[T>G]CG` = "T>G", `CC[T>G]CT` = "T>G", `CC[T>G]GA` = "T>G", `CC[T>G]GC` = "T>G", `CC[T>G]GG` = "T>G", `CC[T>G]GT` = "T>G", `CC[T>G]TA` = "T>G", `CC[T>G]TC` = "T>G", `CC[T>G]TG` = "T>G", `CC[T>G]TT` = "T>G", `CG[C>A]AA` = "C>A", `CG[C>A]AC` = "C>A", `CG[C>A]AG` = "C>A", `CG[C>A]AT` = "C>A", `CG[C>A]CA` = "C>A", `CG[C>A]CC` = "C>A", `CG[C>A]CG` = "C>A", `CG[C>A]CT` = "C>A", `CG[C>A]GA` = "C>A", `CG[C>A]GC` = "C>A", `CG[C>A]GG` = "C>A", `CG[C>A]GT` = "C>A", `CG[C>A]TA` = "C>A", `CG[C>A]TC` = "C>A", `CG[C>A]TG` = "C>A", `CG[C>A]TT` = "C>A", `CG[C>G]AA` = "C>G", `CG[C>G]AC` = "C>G", `CG[C>G]AG` = "C>G", `CG[C>G]AT` = "C>G", `CG[C>G]CA` = "C>G", `CG[C>G]CC` = "C>G", `CG[C>G]CG` = "C>G", `CG[C>G]CT` = "C>G", `CG[C>G]GA` = "C>G", `CG[C>G]GC` = "C>G", `CG[C>G]GG` = "C>G", `CG[C>G]GT` = "C>G", `CG[C>G]TA` = "C>G", `CG[C>G]TC` = "C>G", `CG[C>G]TG` = "C>G", `CG[C>G]TT` = "C>G", `CG[C>T]AA` = "C>T", `CG[C>T]AC` = "C>T", `CG[C>T]AG` = "C>T", `CG[C>T]AT` = "C>T", `CG[C>T]CA` = "C>T", `CG[C>T]CC` = "C>T", `CG[C>T]CG` = "C>T", `CG[C>T]CT` = "C>T", `CG[C>T]GA` = "C>T", `CG[C>T]GC` = "C>T", `CG[C>T]GG` = "C>T", `CG[C>T]GT` = "C>T", `CG[C>T]TA` = "C>T", `CG[C>T]TC` = "C>T", `CG[C>T]TG` = "C>T", `CG[C>T]TT` = "C>T", `CG[T>A]AA` = "T>A", `CG[T>A]AC` = "T>A", `CG[T>A]AG` = "T>A", `CG[T>A]AT` = "T>A", `CG[T>A]CA` = "T>A", `CG[T>A]CC` = "T>A", `CG[T>A]CG` = "T>A", `CG[T>A]CT` = "T>A", `CG[T>A]GA` = "T>A", `CG[T>A]GC` = "T>A", `CG[T>A]GG` = "T>A", `CG[T>A]GT` = "T>A", `CG[T>A]TA` = "T>A", `CG[T>A]TC` = "T>A", `CG[T>A]TG` = "T>A", `CG[T>A]TT` = "T>A", `CG[T>C]AA` = "T>C", `CG[T>C]AC` = "T>C", `CG[T>C]AG` = "T>C", `CG[T>C]AT` = "T>C", `CG[T>C]CA` = "T>C", `CG[T>C]CC` = "T>C", `CG[T>C]CG` = "T>C", `CG[T>C]CT` = "T>C", `CG[T>C]GA` = "T>C", `CG[T>C]GC` = "T>C", `CG[T>C]GG` = "T>C", `CG[T>C]GT` = "T>C", `CG[T>C]TA` = "T>C", `CG[T>C]TC` = "T>C", `CG[T>C]TG` = "T>C", `CG[T>C]TT` = "T>C", `CG[T>G]AA` = "T>G", `CG[T>G]AC` = "T>G", `CG[T>G]AG` = "T>G", `CG[T>G]AT` = "T>G", `CG[T>G]CA` = "T>G", `CG[T>G]CC` = "T>G", `CG[T>G]CG` = "T>G", `CG[T>G]CT` = "T>G", `CG[T>G]GA` = "T>G", `CG[T>G]GC` = "T>G", `CG[T>G]GG` = "T>G", `CG[T>G]GT` = "T>G", `CG[T>G]TA` = "T>G", `CG[T>G]TC` = "T>G", `CG[T>G]TG` = "T>G", `CG[T>G]TT` = "T>G", `CT[C>A]AA` = "C>A", `CT[C>A]AC` = "C>A", `CT[C>A]AG` = "C>A", `CT[C>A]AT` = "C>A", `CT[C>A]CA` = "C>A", `CT[C>A]CC` = "C>A", `CT[C>A]CG` = "C>A", `CT[C>A]CT` = "C>A", `CT[C>A]GA` = "C>A", `CT[C>A]GC` = "C>A", `CT[C>A]GG` = "C>A", `CT[C>A]GT` = "C>A", `CT[C>A]TA` = "C>A", `CT[C>A]TC` = "C>A", `CT[C>A]TG` = "C>A", `CT[C>A]TT` = "C>A", `CT[C>G]AA` = "C>G", `CT[C>G]AC` = "C>G", `CT[C>G]AG` = "C>G", `CT[C>G]AT` = "C>G", `CT[C>G]CA` = "C>G", `CT[C>G]CC` = "C>G", `CT[C>G]CG` = "C>G", `CT[C>G]CT` = "C>G", `CT[C>G]GA` = "C>G", `CT[C>G]GC` = "C>G", `CT[C>G]GG` = "C>G", `CT[C>G]GT` = "C>G", `CT[C>G]TA` = "C>G", `CT[C>G]TC` = "C>G", `CT[C>G]TG` = "C>G", `CT[C>G]TT` = "C>G", `CT[C>T]AA` = "C>T", `CT[C>T]AC` = "C>T", `CT[C>T]AG` = "C>T", `CT[C>T]AT` = "C>T", `CT[C>T]CA` = "C>T", `CT[C>T]CC` = "C>T", `CT[C>T]CG` = "C>T", `CT[C>T]CT` = "C>T", `CT[C>T]GA` = "C>T", `CT[C>T]GC` = "C>T", `CT[C>T]GG` = "C>T", `CT[C>T]GT` = "C>T", `CT[C>T]TA` = "C>T", `CT[C>T]TC` = "C>T", `CT[C>T]TG` = "C>T", `CT[C>T]TT` = "C>T", `CT[T>A]AA` = "T>A", `CT[T>A]AC` = "T>A", `CT[T>A]AG` = "T>A", `CT[T>A]AT` = "T>A", `CT[T>A]CA` = "T>A", `CT[T>A]CC` = "T>A", `CT[T>A]CG` = "T>A", `CT[T>A]CT` = "T>A", `CT[T>A]GA` = "T>A", `CT[T>A]GC` = "T>A", `CT[T>A]GG` = "T>A", `CT[T>A]GT` = "T>A", `CT[T>A]TA` = "T>A", `CT[T>A]TC` = "T>A", `CT[T>A]TG` = "T>A", `CT[T>A]TT` = "T>A", `CT[T>C]AA` = "T>C", `CT[T>C]AC` = "T>C", `CT[T>C]AG` = "T>C", `CT[T>C]AT` = "T>C", `CT[T>C]CA` = "T>C", `CT[T>C]CC` = "T>C", `CT[T>C]CG` = "T>C", `CT[T>C]CT` = "T>C", `CT[T>C]GA` = "T>C", `CT[T>C]GC` = "T>C", `CT[T>C]GG` = "T>C", `CT[T>C]GT` = "T>C", `CT[T>C]TA` = "T>C", `CT[T>C]TC` = "T>C", `CT[T>C]TG` = "T>C", `CT[T>C]TT` = "T>C", `CT[T>G]AA` = "T>G", `CT[T>G]AC` = "T>G", `CT[T>G]AG` = "T>G", `CT[T>G]AT` = "T>G", `CT[T>G]CA` = "T>G", `CT[T>G]CC` = "T>G", `CT[T>G]CG` = "T>G", `CT[T>G]CT` = "T>G", `CT[T>G]GA` = "T>G", `CT[T>G]GC` = "T>G", `CT[T>G]GG` = "T>G", `CT[T>G]GT` = "T>G", `CT[T>G]TA` = "T>G", `CT[T>G]TC` = "T>G", `CT[T>G]TG` = "T>G", `CT[T>G]TT` = "T>G", `GA[C>A]AA` = "C>A", `GA[C>A]AC` = "C>A", `GA[C>A]AG` = "C>A", `GA[C>A]AT` = "C>A", `GA[C>A]CA` = "C>A", `GA[C>A]CC` = "C>A", `GA[C>A]CG` = "C>A", `GA[C>A]CT` = "C>A", `GA[C>A]GA` = "C>A", `GA[C>A]GC` = "C>A", `GA[C>A]GG` = "C>A", `GA[C>A]GT` = "C>A", `GA[C>A]TA` = "C>A", `GA[C>A]TC` = "C>A", `GA[C>A]TG` = "C>A", `GA[C>A]TT` = "C>A", `GA[C>G]AA` = "C>G", `GA[C>G]AC` = "C>G", `GA[C>G]AG` = "C>G", `GA[C>G]AT` = "C>G", `GA[C>G]CA` = "C>G", `GA[C>G]CC` = "C>G", `GA[C>G]CG` = "C>G", `GA[C>G]CT` = "C>G", `GA[C>G]GA` = "C>G", `GA[C>G]GC` = "C>G", `GA[C>G]GG` = "C>G", `GA[C>G]GT` = "C>G", `GA[C>G]TA` = "C>G", `GA[C>G]TC` = "C>G", `GA[C>G]TG` = "C>G", `GA[C>G]TT` = "C>G", `GA[C>T]AA` = "C>T", `GA[C>T]AC` = "C>T", `GA[C>T]AG` = "C>T", `GA[C>T]AT` = "C>T", `GA[C>T]CA` = "C>T", `GA[C>T]CC` = "C>T", `GA[C>T]CG` = "C>T", `GA[C>T]CT` = "C>T", `GA[C>T]GA` = "C>T", `GA[C>T]GC` = "C>T", `GA[C>T]GG` = "C>T", `GA[C>T]GT` = "C>T", `GA[C>T]TA` = "C>T", `GA[C>T]TC` = "C>T", `GA[C>T]TG` = "C>T", `GA[C>T]TT` = "C>T", `GA[T>A]AA` = "T>A", `GA[T>A]AC` = "T>A", `GA[T>A]AG` = "T>A", `GA[T>A]AT` = "T>A", `GA[T>A]CA` = "T>A", `GA[T>A]CC` = "T>A", `GA[T>A]CG` = "T>A", `GA[T>A]CT` = "T>A", `GA[T>A]GA` = "T>A", `GA[T>A]GC` = "T>A", `GA[T>A]GG` = "T>A", `GA[T>A]GT` = "T>A", `GA[T>A]TA` = "T>A", `GA[T>A]TC` = "T>A", `GA[T>A]TG` = "T>A", `GA[T>A]TT` = "T>A", `GA[T>C]AA` = "T>C", `GA[T>C]AC` = "T>C", `GA[T>C]AG` = "T>C", `GA[T>C]AT` = "T>C", `GA[T>C]CA` = "T>C", `GA[T>C]CC` = "T>C", `GA[T>C]CG` = "T>C", `GA[T>C]CT` = "T>C", `GA[T>C]GA` = "T>C", `GA[T>C]GC` = "T>C", `GA[T>C]GG` = "T>C", `GA[T>C]GT` = "T>C", `GA[T>C]TA` = "T>C", `GA[T>C]TC` = "T>C", `GA[T>C]TG` = "T>C", `GA[T>C]TT` = "T>C", `GA[T>G]AA` = "T>G", `GA[T>G]AC` = "T>G", `GA[T>G]AG` = "T>G", `GA[T>G]AT` = "T>G", `GA[T>G]CA` = "T>G", `GA[T>G]CC` = "T>G", `GA[T>G]CG` = "T>G", `GA[T>G]CT` = "T>G", `GA[T>G]GA` = "T>G", `GA[T>G]GC` = "T>G", `GA[T>G]GG` = "T>G", `GA[T>G]GT` = "T>G", `GA[T>G]TA` = "T>G", `GA[T>G]TC` = "T>G", `GA[T>G]TG` = "T>G", `GA[T>G]TT` = "T>G", `GC[C>A]AA` = "C>A", `GC[C>A]AC` = "C>A", `GC[C>A]AG` = "C>A", `GC[C>A]AT` = "C>A", `GC[C>A]CA` = "C>A", `GC[C>A]CC` = "C>A", `GC[C>A]CG` = "C>A", `GC[C>A]CT` = "C>A", `GC[C>A]GA` = "C>A", `GC[C>A]GC` = "C>A", `GC[C>A]GG` = "C>A", `GC[C>A]GT` = "C>A", `GC[C>A]TA` = "C>A", `GC[C>A]TC` = "C>A", `GC[C>A]TG` = "C>A", `GC[C>A]TT` = "C>A", `GC[C>G]AA` = "C>G", `GC[C>G]AC` = "C>G", `GC[C>G]AG` = "C>G", `GC[C>G]AT` = "C>G", `GC[C>G]CA` = "C>G", `GC[C>G]CC` = "C>G", `GC[C>G]CG` = "C>G", `GC[C>G]CT` = "C>G", `GC[C>G]GA` = "C>G", `GC[C>G]GC` = "C>G", `GC[C>G]GG` = "C>G", `GC[C>G]GT` = "C>G", `GC[C>G]TA` = "C>G", `GC[C>G]TC` = "C>G", `GC[C>G]TG` = "C>G", `GC[C>G]TT` = "C>G", `GC[C>T]AA` = "C>T", `GC[C>T]AC` = "C>T", `GC[C>T]AG` = "C>T", `GC[C>T]AT` = "C>T", `GC[C>T]CA` = "C>T", `GC[C>T]CC` = "C>T", `GC[C>T]CG` = "C>T", `GC[C>T]CT` = "C>T", `GC[C>T]GA` = "C>T", `GC[C>T]GC` = "C>T", `GC[C>T]GG` = "C>T", `GC[C>T]GT` = "C>T", `GC[C>T]TA` = "C>T", `GC[C>T]TC` = "C>T", `GC[C>T]TG` = "C>T", `GC[C>T]TT` = "C>T", `GC[T>A]AA` = "T>A", `GC[T>A]AC` = "T>A", `GC[T>A]AG` = "T>A", `GC[T>A]AT` = "T>A", `GC[T>A]CA` = "T>A", `GC[T>A]CC` = "T>A", `GC[T>A]CG` = "T>A", `GC[T>A]CT` = "T>A", `GC[T>A]GA` = "T>A", `GC[T>A]GC` = "T>A", `GC[T>A]GG` = "T>A", `GC[T>A]GT` = "T>A", `GC[T>A]TA` = "T>A", `GC[T>A]TC` = "T>A", `GC[T>A]TG` = "T>A", `GC[T>A]TT` = "T>A", `GC[T>C]AA` = "T>C", `GC[T>C]AC` = "T>C", `GC[T>C]AG` = "T>C", `GC[T>C]AT` = "T>C", `GC[T>C]CA` = "T>C", `GC[T>C]CC` = "T>C", `GC[T>C]CG` = "T>C", `GC[T>C]CT` = "T>C", `GC[T>C]GA` = "T>C", `GC[T>C]GC` = "T>C", `GC[T>C]GG` = "T>C", `GC[T>C]GT` = "T>C", `GC[T>C]TA` = "T>C", `GC[T>C]TC` = "T>C", `GC[T>C]TG` = "T>C", `GC[T>C]TT` = "T>C", `GC[T>G]AA` = "T>G", `GC[T>G]AC` = "T>G", `GC[T>G]AG` = "T>G", `GC[T>G]AT` = "T>G", `GC[T>G]CA` = "T>G", `GC[T>G]CC` = "T>G", `GC[T>G]CG` = "T>G", `GC[T>G]CT` = "T>G", `GC[T>G]GA` = "T>G", `GC[T>G]GC` = "T>G", `GC[T>G]GG` = "T>G", `GC[T>G]GT` = "T>G", `GC[T>G]TA` = "T>G", `GC[T>G]TC` = "T>G", `GC[T>G]TG` = "T>G", `GC[T>G]TT` = "T>G", `GG[C>A]AA` = "C>A", `GG[C>A]AC` = "C>A", `GG[C>A]AG` = "C>A", `GG[C>A]AT` = "C>A", `GG[C>A]CA` = "C>A", `GG[C>A]CC` = "C>A", `GG[C>A]CG` = "C>A", `GG[C>A]CT` = "C>A", `GG[C>A]GA` = "C>A", `GG[C>A]GC` = "C>A", `GG[C>A]GG` = "C>A", `GG[C>A]GT` = "C>A", `GG[C>A]TA` = "C>A", `GG[C>A]TC` = "C>A", `GG[C>A]TG` = "C>A", `GG[C>A]TT` = "C>A", `GG[C>G]AA` = "C>G", `GG[C>G]AC` = "C>G", `GG[C>G]AG` = "C>G", `GG[C>G]AT` = "C>G", `GG[C>G]CA` = "C>G", `GG[C>G]CC` = "C>G", `GG[C>G]CG` = "C>G", `GG[C>G]CT` = "C>G", `GG[C>G]GA` = "C>G", `GG[C>G]GC` = "C>G", `GG[C>G]GG` = "C>G", `GG[C>G]GT` = "C>G", `GG[C>G]TA` = "C>G", `GG[C>G]TC` = "C>G", `GG[C>G]TG` = "C>G", `GG[C>G]TT` = "C>G", `GG[C>T]AA` = "C>T", `GG[C>T]AC` = "C>T", `GG[C>T]AG` = "C>T", `GG[C>T]AT` = "C>T", `GG[C>T]CA` = "C>T", `GG[C>T]CC` = "C>T", `GG[C>T]CG` = "C>T", `GG[C>T]CT` = "C>T", `GG[C>T]GA` = "C>T", `GG[C>T]GC` = "C>T", `GG[C>T]GG` = "C>T", `GG[C>T]GT` = "C>T", `GG[C>T]TA` = "C>T", `GG[C>T]TC` = "C>T", `GG[C>T]TG` = "C>T", `GG[C>T]TT` = "C>T", `GG[T>A]AA` = "T>A", `GG[T>A]AC` = "T>A", `GG[T>A]AG` = "T>A", `GG[T>A]AT` = "T>A", `GG[T>A]CA` = "T>A", `GG[T>A]CC` = "T>A", `GG[T>A]CG` = "T>A", `GG[T>A]CT` = "T>A", `GG[T>A]GA` = "T>A", `GG[T>A]GC` = "T>A", `GG[T>A]GG` = "T>A", `GG[T>A]GT` = "T>A", `GG[T>A]TA` = "T>A", `GG[T>A]TC` = "T>A", `GG[T>A]TG` = "T>A", `GG[T>A]TT` = "T>A", `GG[T>C]AA` = "T>C", `GG[T>C]AC` = "T>C", `GG[T>C]AG` = "T>C", `GG[T>C]AT` = "T>C", `GG[T>C]CA` = "T>C", `GG[T>C]CC` = "T>C", `GG[T>C]CG` = "T>C", `GG[T>C]CT` = "T>C", `GG[T>C]GA` = "T>C", `GG[T>C]GC` = "T>C", `GG[T>C]GG` = "T>C", `GG[T>C]GT` = "T>C", `GG[T>C]TA` = "T>C", `GG[T>C]TC` = "T>C", `GG[T>C]TG` = "T>C", `GG[T>C]TT` = "T>C", `GG[T>G]AA` = "T>G", `GG[T>G]AC` = "T>G", `GG[T>G]AG` = "T>G", `GG[T>G]AT` = "T>G", `GG[T>G]CA` = "T>G", `GG[T>G]CC` = "T>G", `GG[T>G]CG` = "T>G", `GG[T>G]CT` = "T>G", `GG[T>G]GA` = "T>G", `GG[T>G]GC` = "T>G", `GG[T>G]GG` = "T>G", `GG[T>G]GT` = "T>G", `GG[T>G]TA` = "T>G", `GG[T>G]TC` = "T>G", `GG[T>G]TG` = "T>G", `GG[T>G]TT` = "T>G", `GT[C>A]AA` = "C>A", `GT[C>A]AC` = "C>A", `GT[C>A]AG` = "C>A", `GT[C>A]AT` = "C>A", `GT[C>A]CA` = "C>A", `GT[C>A]CC` = "C>A", `GT[C>A]CG` = "C>A", `GT[C>A]CT` = "C>A", `GT[C>A]GA` = "C>A", `GT[C>A]GC` = "C>A", `GT[C>A]GG` = "C>A", `GT[C>A]GT` = "C>A", `GT[C>A]TA` = "C>A", `GT[C>A]TC` = "C>A", `GT[C>A]TG` = "C>A", `GT[C>A]TT` = "C>A", `GT[C>G]AA` = "C>G", `GT[C>G]AC` = "C>G", `GT[C>G]AG` = "C>G", `GT[C>G]AT` = "C>G", `GT[C>G]CA` = "C>G", `GT[C>G]CC` = "C>G", `GT[C>G]CG` = "C>G", `GT[C>G]CT` = "C>G", `GT[C>G]GA` = "C>G", `GT[C>G]GC` = "C>G", `GT[C>G]GG` = "C>G", `GT[C>G]GT` = "C>G", `GT[C>G]TA` = "C>G", `GT[C>G]TC` = "C>G", `GT[C>G]TG` = "C>G", `GT[C>G]TT` = "C>G", `GT[C>T]AA` = "C>T", `GT[C>T]AC` = "C>T", `GT[C>T]AG` = "C>T", `GT[C>T]AT` = "C>T", `GT[C>T]CA` = "C>T", `GT[C>T]CC` = "C>T", `GT[C>T]CG` = "C>T", `GT[C>T]CT` = "C>T", `GT[C>T]GA` = "C>T", `GT[C>T]GC` = "C>T", `GT[C>T]GG` = "C>T", `GT[C>T]GT` = "C>T", `GT[C>T]TA` = "C>T", `GT[C>T]TC` = "C>T", `GT[C>T]TG` = "C>T", `GT[C>T]TT` = "C>T", `GT[T>A]AA` = "T>A", `GT[T>A]AC` = "T>A", `GT[T>A]AG` = "T>A", `GT[T>A]AT` = "T>A", `GT[T>A]CA` = "T>A", `GT[T>A]CC` = "T>A", `GT[T>A]CG` = "T>A", `GT[T>A]CT` = "T>A", `GT[T>A]GA` = "T>A", `GT[T>A]GC` = "T>A", `GT[T>A]GG` = "T>A", `GT[T>A]GT` = "T>A", `GT[T>A]TA` = "T>A", `GT[T>A]TC` = "T>A", `GT[T>A]TG` = "T>A", `GT[T>A]TT` = "T>A", `GT[T>C]AA` = "T>C", `GT[T>C]AC` = "T>C", `GT[T>C]AG` = "T>C", `GT[T>C]AT` = "T>C", `GT[T>C]CA` = "T>C", `GT[T>C]CC` = "T>C", `GT[T>C]CG` = "T>C", `GT[T>C]CT` = "T>C", `GT[T>C]GA` = "T>C", `GT[T>C]GC` = "T>C", `GT[T>C]GG` = "T>C", `GT[T>C]GT` = "T>C", `GT[T>C]TA` = "T>C", `GT[T>C]TC` = "T>C", `GT[T>C]TG` = "T>C", `GT[T>C]TT` = "T>C", `GT[T>G]AA` = "T>G", `GT[T>G]AC` = "T>G", `GT[T>G]AG` = "T>G", `GT[T>G]AT` = "T>G", `GT[T>G]CA` = "T>G", `GT[T>G]CC` = "T>G", `GT[T>G]CG` = "T>G", `GT[T>G]CT` = "T>G", `GT[T>G]GA` = "T>G", `GT[T>G]GC` = "T>G", `GT[T>G]GG` = "T>G", `GT[T>G]GT` = "T>G", `GT[T>G]TA` = "T>G", `GT[T>G]TC` = "T>G", `GT[T>G]TG` = "T>G", `GT[T>G]TT` = "T>G", `TA[C>A]AA` = "C>A", `TA[C>A]AC` = "C>A", `TA[C>A]AG` = "C>A", `TA[C>A]AT` = "C>A", `TA[C>A]CA` = "C>A", `TA[C>A]CC` = "C>A", `TA[C>A]CG` = "C>A", `TA[C>A]CT` = "C>A", `TA[C>A]GA` = "C>A", `TA[C>A]GC` = "C>A", `TA[C>A]GG` = "C>A", `TA[C>A]GT` = "C>A", `TA[C>A]TA` = "C>A", `TA[C>A]TC` = "C>A", `TA[C>A]TG` = "C>A", `TA[C>A]TT` = "C>A", `TA[C>G]AA` = "C>G", `TA[C>G]AC` = "C>G", `TA[C>G]AG` = "C>G", `TA[C>G]AT` = "C>G", `TA[C>G]CA` = "C>G", `TA[C>G]CC` = "C>G", `TA[C>G]CG` = "C>G", `TA[C>G]CT` = "C>G", `TA[C>G]GA` = "C>G", `TA[C>G]GC` = "C>G", `TA[C>G]GG` = "C>G", `TA[C>G]GT` = "C>G", `TA[C>G]TA` = "C>G", `TA[C>G]TC` = "C>G", `TA[C>G]TG` = "C>G", `TA[C>G]TT` = "C>G", `TA[C>T]AA` = "C>T", `TA[C>T]AC` = "C>T", `TA[C>T]AG` = "C>T", `TA[C>T]AT` = "C>T", `TA[C>T]CA` = "C>T", `TA[C>T]CC` = "C>T", `TA[C>T]CG` = "C>T", `TA[C>T]CT` = "C>T", `TA[C>T]GA` = "C>T", `TA[C>T]GC` = "C>T", `TA[C>T]GG` = "C>T", `TA[C>T]GT` = "C>T", `TA[C>T]TA` = "C>T", `TA[C>T]TC` = "C>T", `TA[C>T]TG` = "C>T", `TA[C>T]TT` = "C>T", `TA[T>A]AA` = "T>A", `TA[T>A]AC` = "T>A", `TA[T>A]AG` = "T>A", `TA[T>A]AT` = "T>A", `TA[T>A]CA` = "T>A", `TA[T>A]CC` = "T>A", `TA[T>A]CG` = "T>A", `TA[T>A]CT` = "T>A", `TA[T>A]GA` = "T>A", `TA[T>A]GC` = "T>A", `TA[T>A]GG` = "T>A", `TA[T>A]GT` = "T>A", `TA[T>A]TA` = "T>A", `TA[T>A]TC` = "T>A", `TA[T>A]TG` = "T>A", `TA[T>A]TT` = "T>A", `TA[T>C]AA` = "T>C", `TA[T>C]AC` = "T>C", `TA[T>C]AG` = "T>C", `TA[T>C]AT` = "T>C", `TA[T>C]CA` = "T>C", `TA[T>C]CC` = "T>C", `TA[T>C]CG` = "T>C", `TA[T>C]CT` = "T>C", `TA[T>C]GA` = "T>C", `TA[T>C]GC` = "T>C", `TA[T>C]GG` = "T>C", `TA[T>C]GT` = "T>C", `TA[T>C]TA` = "T>C", `TA[T>C]TC` = "T>C", `TA[T>C]TG` = "T>C", `TA[T>C]TT` = "T>C", `TA[T>G]AA` = "T>G", `TA[T>G]AC` = "T>G", `TA[T>G]AG` = "T>G", `TA[T>G]AT` = "T>G", `TA[T>G]CA` = "T>G", `TA[T>G]CC` = "T>G", `TA[T>G]CG` = "T>G", `TA[T>G]CT` = "T>G", `TA[T>G]GA` = "T>G", `TA[T>G]GC` = "T>G", `TA[T>G]GG` = "T>G", `TA[T>G]GT` = "T>G", `TA[T>G]TA` = "T>G", `TA[T>G]TC` = "T>G", `TA[T>G]TG` = "T>G", `TA[T>G]TT` = "T>G", `TC[C>A]AA` = "C>A", `TC[C>A]AC` = "C>A", `TC[C>A]AG` = "C>A", `TC[C>A]AT` = "C>A", `TC[C>A]CA` = "C>A", `TC[C>A]CC` = "C>A", `TC[C>A]CG` = "C>A", `TC[C>A]CT` = "C>A", `TC[C>A]GA` = "C>A", `TC[C>A]GC` = "C>A", `TC[C>A]GG` = "C>A", `TC[C>A]GT` = "C>A", `TC[C>A]TA` = "C>A", `TC[C>A]TC` = "C>A", `TC[C>A]TG` = "C>A", `TC[C>A]TT` = "C>A", `TC[C>G]AA` = "C>G", `TC[C>G]AC` = "C>G", `TC[C>G]AG` = "C>G", `TC[C>G]AT` = "C>G", `TC[C>G]CA` = "C>G", `TC[C>G]CC` = "C>G", `TC[C>G]CG` = "C>G", `TC[C>G]CT` = "C>G", `TC[C>G]GA` = "C>G", `TC[C>G]GC` = "C>G", `TC[C>G]GG` = "C>G", `TC[C>G]GT` = "C>G", `TC[C>G]TA` = "C>G", `TC[C>G]TC` = "C>G", `TC[C>G]TG` = "C>G", `TC[C>G]TT` = "C>G", `TC[C>T]AA` = "C>T", `TC[C>T]AC` = "C>T", `TC[C>T]AG` = "C>T", `TC[C>T]AT` = "C>T", `TC[C>T]CA` = "C>T", `TC[C>T]CC` = "C>T", `TC[C>T]CG` = "C>T", `TC[C>T]CT` = "C>T", `TC[C>T]GA` = "C>T", `TC[C>T]GC` = "C>T", `TC[C>T]GG` = "C>T", `TC[C>T]GT` = "C>T", `TC[C>T]TA` = "C>T", `TC[C>T]TC` = "C>T", `TC[C>T]TG` = "C>T", `TC[C>T]TT` = "C>T", `TC[T>A]AA` = "T>A", `TC[T>A]AC` = "T>A", `TC[T>A]AG` = "T>A", `TC[T>A]AT` = "T>A", `TC[T>A]CA` = "T>A", `TC[T>A]CC` = "T>A", `TC[T>A]CG` = "T>A", `TC[T>A]CT` = "T>A", `TC[T>A]GA` = "T>A", `TC[T>A]GC` = "T>A", `TC[T>A]GG` = "T>A", `TC[T>A]GT` = "T>A", `TC[T>A]TA` = "T>A", `TC[T>A]TC` = "T>A", `TC[T>A]TG` = "T>A", `TC[T>A]TT` = "T>A", `TC[T>C]AA` = "T>C", `TC[T>C]AC` = "T>C", `TC[T>C]AG` = "T>C", `TC[T>C]AT` = "T>C", `TC[T>C]CA` = "T>C", `TC[T>C]CC` = "T>C", `TC[T>C]CG` = "T>C", `TC[T>C]CT` = "T>C", `TC[T>C]GA` = "T>C", `TC[T>C]GC` = "T>C", `TC[T>C]GG` = "T>C", `TC[T>C]GT` = "T>C", `TC[T>C]TA` = "T>C", `TC[T>C]TC` = "T>C", `TC[T>C]TG` = "T>C", `TC[T>C]TT` = "T>C", `TC[T>G]AA` = "T>G", `TC[T>G]AC` = "T>G", `TC[T>G]AG` = "T>G", `TC[T>G]AT` = "T>G", `TC[T>G]CA` = "T>G", `TC[T>G]CC` = "T>G", `TC[T>G]CG` = "T>G", `TC[T>G]CT` = "T>G", `TC[T>G]GA` = "T>G", `TC[T>G]GC` = "T>G", `TC[T>G]GG` = "T>G", `TC[T>G]GT` = "T>G", `TC[T>G]TA` = "T>G", `TC[T>G]TC` = "T>G", `TC[T>G]TG` = "T>G", `TC[T>G]TT` = "T>G", `TG[C>A]AA` = "C>A", `TG[C>A]AC` = "C>A", `TG[C>A]AG` = "C>A", `TG[C>A]AT` = "C>A", `TG[C>A]CA` = "C>A", `TG[C>A]CC` = "C>A", `TG[C>A]CG` = "C>A", `TG[C>A]CT` = "C>A", `TG[C>A]GA` = "C>A", `TG[C>A]GC` = "C>A", `TG[C>A]GG` = "C>A", `TG[C>A]GT` = "C>A", `TG[C>A]TA` = "C>A", `TG[C>A]TC` = "C>A", `TG[C>A]TG` = "C>A", `TG[C>A]TT` = "C>A", `TG[C>G]AA` = "C>G", `TG[C>G]AC` = "C>G", `TG[C>G]AG` = "C>G", `TG[C>G]AT` = "C>G", `TG[C>G]CA` = "C>G", `TG[C>G]CC` = "C>G", `TG[C>G]CG` = "C>G", `TG[C>G]CT` = "C>G", `TG[C>G]GA` = "C>G", `TG[C>G]GC` = "C>G", `TG[C>G]GG` = "C>G", `TG[C>G]GT` = "C>G", `TG[C>G]TA` = "C>G", `TG[C>G]TC` = "C>G", `TG[C>G]TG` = "C>G", `TG[C>G]TT` = "C>G", `TG[C>T]AA` = "C>T", `TG[C>T]AC` = "C>T", `TG[C>T]AG` = "C>T", `TG[C>T]AT` = "C>T", `TG[C>T]CA` = "C>T", `TG[C>T]CC` = "C>T", `TG[C>T]CG` = "C>T", `TG[C>T]CT` = "C>T", `TG[C>T]GA` = "C>T", `TG[C>T]GC` = "C>T", `TG[C>T]GG` = "C>T", `TG[C>T]GT` = "C>T", `TG[C>T]TA` = "C>T", `TG[C>T]TC` = "C>T", `TG[C>T]TG` = "C>T", `TG[C>T]TT` = "C>T", `TG[T>A]AA` = "T>A", `TG[T>A]AC` = "T>A", `TG[T>A]AG` = "T>A", `TG[T>A]AT` = "T>A", `TG[T>A]CA` = "T>A", `TG[T>A]CC` = "T>A", `TG[T>A]CG` = "T>A", `TG[T>A]CT` = "T>A", `TG[T>A]GA` = "T>A", `TG[T>A]GC` = "T>A", `TG[T>A]GG` = "T>A", `TG[T>A]GT` = "T>A", `TG[T>A]TA` = "T>A", `TG[T>A]TC` = "T>A", `TG[T>A]TG` = "T>A", `TG[T>A]TT` = "T>A", `TG[T>C]AA` = "T>C", `TG[T>C]AC` = "T>C", `TG[T>C]AG` = "T>C", `TG[T>C]AT` = "T>C", `TG[T>C]CA` = "T>C", `TG[T>C]CC` = "T>C", `TG[T>C]CG` = "T>C", `TG[T>C]CT` = "T>C", `TG[T>C]GA` = "T>C", `TG[T>C]GC` = "T>C", `TG[T>C]GG` = "T>C", `TG[T>C]GT` = "T>C", `TG[T>C]TA` = "T>C", `TG[T>C]TC` = "T>C", `TG[T>C]TG` = "T>C", `TG[T>C]TT` = "T>C", `TG[T>G]AA` = "T>G", `TG[T>G]AC` = "T>G", `TG[T>G]AG` = "T>G", `TG[T>G]AT` = "T>G", `TG[T>G]CA` = "T>G", `TG[T>G]CC` = "T>G", `TG[T>G]CG` = "T>G", `TG[T>G]CT` = "T>G", `TG[T>G]GA` = "T>G", `TG[T>G]GC` = "T>G", `TG[T>G]GG` = "T>G", `TG[T>G]GT` = "T>G", `TG[T>G]TA` = "T>G", `TG[T>G]TC` = "T>G", `TG[T>G]TG` = "T>G", `TG[T>G]TT` = "T>G", `TT[C>A]AA` = "C>A", `TT[C>A]AC` = "C>A", `TT[C>A]AG` = "C>A", `TT[C>A]AT` = "C>A", `TT[C>A]CA` = "C>A", `TT[C>A]CC` = "C>A", `TT[C>A]CG` = "C>A", `TT[C>A]CT` = "C>A", `TT[C>A]GA` = "C>A", `TT[C>A]GC` = "C>A", `TT[C>A]GG` = "C>A", `TT[C>A]GT` = "C>A", `TT[C>A]TA` = "C>A", `TT[C>A]TC` = "C>A", `TT[C>A]TG` = "C>A", `TT[C>A]TT` = "C>A", `TT[C>G]AA` = "C>G", `TT[C>G]AC` = "C>G", `TT[C>G]AG` = "C>G", `TT[C>G]AT` = "C>G", `TT[C>G]CA` = "C>G", `TT[C>G]CC` = "C>G", `TT[C>G]CG` = "C>G", `TT[C>G]CT` = "C>G", `TT[C>G]GA` = "C>G", `TT[C>G]GC` = "C>G", `TT[C>G]GG` = "C>G", `TT[C>G]GT` = "C>G", `TT[C>G]TA` = "C>G", `TT[C>G]TC` = "C>G", `TT[C>G]TG` = "C>G", `TT[C>G]TT` = "C>G", `TT[C>T]AA` = "C>T", `TT[C>T]AC` = "C>T", `TT[C>T]AG` = "C>T", `TT[C>T]AT` = "C>T", `TT[C>T]CA` = "C>T", `TT[C>T]CC` = "C>T", `TT[C>T]CG` = "C>T", `TT[C>T]CT` = "C>T", `TT[C>T]GA` = "C>T", `TT[C>T]GC` = "C>T", `TT[C>T]GG` = "C>T", `TT[C>T]GT` = "C>T", `TT[C>T]TA` = "C>T", `TT[C>T]TC` = "C>T", `TT[C>T]TG` = "C>T", `TT[C>T]TT` = "C>T", `TT[T>A]AA` = "T>A", `TT[T>A]AC` = "T>A", `TT[T>A]AG` = "T>A", `TT[T>A]AT` = "T>A", `TT[T>A]CA` = "T>A", `TT[T>A]CC` = "T>A", `TT[T>A]CG` = "T>A", `TT[T>A]CT` = "T>A", `TT[T>A]GA` = "T>A", `TT[T>A]GC` = "T>A", `TT[T>A]GG` = "T>A", `TT[T>A]GT` = "T>A", `TT[T>A]TA` = "T>A", `TT[T>A]TC` = "T>A", `TT[T>A]TG` = "T>A", `TT[T>A]TT` = "T>A", `TT[T>C]AA` = "T>C", `TT[T>C]AC` = "T>C", `TT[T>C]AG` = "T>C", `TT[T>C]AT` = "T>C", `TT[T>C]CA` = "T>C", `TT[T>C]CC` = "T>C", `TT[T>C]CG` = "T>C", `TT[T>C]CT` = "T>C", `TT[T>C]GA` = "T>C", `TT[T>C]GC` = "T>C", `TT[T>C]GG` = "T>C", `TT[T>C]GT` = "T>C", `TT[T>C]TA` = "T>C", `TT[T>C]TC` = "T>C", `TT[T>C]TG` = "T>C", `TT[T>C]TT` = "T>C", `TT[T>G]AA` = "T>G", `TT[T>G]AC` = "T>G", `TT[T>G]AG` = "T>G", `TT[T>G]AT` = "T>G", `TT[T>G]CA` = "T>G", `TT[T>G]CC` = "T>G", `TT[T>G]CG` = "T>G", `TT[T>G]CT` = "T>G", `TT[T>G]GA` = "T>G", `TT[T>G]GC` = "T>G", `TT[T>G]GG` = "T>G", `TT[T>G]GT` = "T>G", `TT[T>G]TA` = "T>G", `TT[T>G]TC` = "T>G", `TT[T>G]TG` = "T>G", `TT[T>G]TT` = "T>G")
}

cosmic_dbs78_channel_to_type <- function() {
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

cosmic_id83_channel_to_type <- function() {
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

cosmic_sv32_channel_to_type <- function() {
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

cosmic_rna_sbs192_channel_to_type <- function() {
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

cosmic_sv38_channel_to_type <- function() {
  c(`clustered:del:<1Kb` = "clustered:del", `clustered:del:>10Mb` = "clustered:del", `clustered:del:1-10Kb` = "clustered:del", `clustered:del:10-100Kb` = "clustered:del", `clustered:del:100Kb-1Mb` = "clustered:del", `clustered:del:1Mb-10Mb` = "clustered:del", `clustered:inv:<1Kb` = "clustered:inv", `clustered:inv:>10Mb` = "clustered:inv", `clustered:inv:1-10Kb` = "clustered:inv", `clustered:inv:10-100Kb` = "clustered:inv", `clustered:inv:100Kb-1Mb` = "clustered:inv", `clustered:inv:1Mb-10Mb` = "clustered:inv", `clustered:tds:<1Kb` = "clustered:tds", `clustered:tds:>10Mb` = "clustered:tds", `clustered:tds:1-10Kb` = "clustered:tds", `clustered:tds:10-100Kb` = "clustered:tds", `clustered:tds:100Kb-1Mb` = "clustered:tds", `clustered:tds:1Mb-10Mb` = "clustered:tds", `clustered:trans` = "clustered:trans", `non-clustered:del:<1Kb` = "non-clustered:del", `non-clustered:del:>10Mb` = "non-clustered:del", `non-clustered:del:1-10Kb` = "non-clustered:del", `non-clustered:del:10-100Kb` = "non-clustered:del", `non-clustered:del:100Kb-1Mb` = "non-clustered:del", `non-clustered:del:1Mb-10Mb` = "non-clustered:del", `non-clustered:inv:<1Kb` = "non-clustered:inv", `non-clustered:inv:>10Mb` = "non-clustered:inv", `non-clustered:inv:1-10Kb` = "non-clustered:inv", `non-clustered:inv:10-100Kb` = "non-clustered:inv", `non-clustered:inv:100Kb-1Mb` = "non-clustered:inv", `non-clustered:inv:1Mb-10Mb` = "non-clustered:inv", `non-clustered:tds:<1Kb` = "non-clustered:tds", `non-clustered:tds:>10Mb` = "non-clustered:tds", `non-clustered:tds:1-10Kb` = "non-clustered:tds", `non-clustered:tds:10-100Kb` = "non-clustered:tds", `non-clustered:tds:100Kb-1Mb` = "non-clustered:tds", `non-clustered:tds:1Mb-10Mb` = "non-clustered:tds", `non-clustered:trans` = "non-clustered:trans")
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

  # Convert Sigstash (cosmic-style) Channel Names to Sigminer
  first_sig_channel_order <- sig_convert_channel_name(first_sig_channel_order, from = "cosmic", to = "sigminer")

  rownames(df_wide) <- first_sig_channel_order


  return(df_wide)
}



# Channel Name Conversions ------------------------------------------------
#' Convert Channel Names Between 'cosmic' and 'sigminer' Formats
#'
#' This function converts channel names between 'cosmic' and 'sigminer' formats. It modifies the capitalization of 'Kb' and replaces certain characters in the channel names as per the required format.
#'
#' @param channel A character vector of channel names to be converted.
#' @param from A character string indicating the format to convert from. Should be either "cosmic" or "sigminer". Defaults to "cosmic".
#' @param to A character string indicating the format to convert to. Should be either "sigminer" or "cosmic". Defaults to "sigminer".
#'
#' @return A character vector with the converted channel names.
#' @export
#'
#' @examples
#' # Convert from 'cosmic' to 'sigminer'
#' sig_convert_channel_name("clustered_del_1-10Kb", from = "cosmic", to = "sigminer")
#'
#' # Convert from 'sigminer' to 'cosmic'
#' sig_convert_channel_name("clustered:del:1-10Kb", from = "sigminer", to = "cosmic")
#'
sig_convert_channel_name <- function(channel, from = c("cosmic", "sigminer"), to = c("sigminer", "cosmic")) {
  from <- rlang::arg_match(from)
  to <- rlang::arg_match(to)

  if (from == "cosmic" & to == "sigminer") {
    # Capitalise kb to Kb in copynumber channel names
    channel <- sub(x = channel, "([0-9])kb", "\\1Kb")

    # Convert _ to : in SV channel names
    channel <- gsub(x = channel, "_", ":")
  } else if (from == "sigminer" & to == "cosmic") {
    # Un-capitalise Kb to kb in copynumber channel names
    channel <- sub(x = channel, "([0-9])Kb", "\\1kb")

    # Convert ':' to '_' in SV channel names
    channel <- gsub(x = channel, ":", "_")
  } else {
    stop("Channel name conversion from ", from, " to ", to, " has not yet been implemented. Please create a new issue on the sigstash github page")
  }
}
