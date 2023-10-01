# Convert: COSMIC > SIGVERSE DataFrame
#
# Cosmic Style to sigverse style signature collection files.
# Note you'll still need to manually add a 'type' column convert signatures
#
# @param data cosmic style dataframe.
# Has header line, 1 column for each signature, and first column corresponds to the channel.
# Each Sig column is a numeric vector describing the fraction of mutations in each group
#
# @return list of data.frames
#
# @examples
# orig_file = 'COSMIC_v3.3.1_SBS_GRCh37.txt'
# path = system.file("original_files", orig_file, package = "sigstash")
# data = utils::read.csv(path, header = TRUE, sep = "\t")
#
# # Convert to sigstash-style
# df_sigs = sig_cosmic_to_sigstash(data, sigclass = 'SNV')
#
# # Write output
# outfile = paste0(tools::file_path_sans_ext(basename(orig_file)), '.csv')
# utils::write.csv(df_sigs, file = paste0("inst/reference_signatures/", outfile), row.names = FALSE, quote = TRUE)
sig_cosmic_to_sigstash <- function(data, sigclass = c("SBS", "ID", "CN", "DBS")){
  requireNamespace("rlang", quietly = TRUE)

  sigclass = rlang::arg_match(sigclass)

  signames = colnames(data)[-1]
  channels = data[[1]]
  types = channel2type(channels, sigclass)

  # For each signature ...
  ls_signatures = lapply(X = signames, FUN = \(sig){

    # Create a data.frame with just Type (renamed to channel) and fractions
    df_sig = data[,c('Type', sig)]
    colnames(df_sig) <- c('channel', 'fraction')

    # Add sigverse defined 'type' column, and 'signature' column describing name
    df_sig[['type']] <- types
    df_sig[['signature']] <- sig
    df_sig <- df_sig[, c('signature', 'type', 'channel', 'fraction')]

    # Convert to tibble
    df_sig = tibble::tibble(df_sig)

    # Return 3 column tibble (channel, type, fraction)
    return(df_sig)

    })

  df_signatures = do.call('rbind', ls_signatures)
  return(df_signatures)
  #names(ls_signatures) <- signames

  #return(ls_signatures)
}

channel2type <- function(channel, sigclass = c("SBS", "ID", "CN", "DBS")){
  requireNamespace("rlang", quietly = TRUE)
  sigclass = rlang::arg_match(sigclass)

  if(sigclass == "SBS")
    vec_channel2type <- cosmic_sbs_channel_to_type()
  else if(sigclass == "CN")
    vec_channel2type <- cosmic_cn_channel_to_type()
  else if(sigclass == "ID")
    vec_channel2type <- cosmic_id_channel_to_type()
  else if(sigclass == "DBS")
    vec_channel2type <- cosmic_dbs_channel_to_type()
  else
    stop('unexpected sigclass')

  types <- vec_channel2type[match(channel, names(vec_channel2type))]
  return(types)

}

cosmic_sbs_channel_to_type <- function(){

  c(`A[C>A]A` = "C>A", `A[C>A]C` = "C>A", `A[C>A]G` = "C>A", `A[C>A]T` = "C>A",
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


cosmic_cn_channel_to_type <- function(){
  c(`0:homdel:0-100kb` = "0", `0:homdel:100kb-1Mb` = "0", `0:homdel:>1Mb` = "0",
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

cosmic_dbs_channel_to_type <- function(){
  c(`AC>CA` = "AC>CA", `AC>CG` = "AC>CG", `AC>CT` = "AC>CT", `AC>GA` = "AC>GA",
    `AC>GG` = "AC>GG", `AC>GT` = "AC>GT", `AC>TA` = "AC>TA", `AC>TG` = "AC>TG",
    `AC>TT` = "AC>TT", `AT>CA` = "AT>CA", `AT>CC` = "AT>CC", `AT>CG` = "AT>CG",
    `AT>GA` = "AT>GA", `AT>GC` = "AT>GC", `AT>TA` = "AT>TA", `CC>AA` = "CC>AA",
    `CC>AG` = "CC>AG", `CC>AT` = "CC>AT", `CC>GA` = "CC>GA", `CC>GG` = "CC>GG",
    `CC>GT` = "CC>GT", `CC>TA` = "CC>TA", `CC>TG` = "CC>TG", `CC>TT` = "CC>TT",
    `CG>AT` = "CG>AT", `CG>GC` = "CG>GC", `CG>GT` = "CG>GT", `CG>TA` = "CG>TA",
    `CG>TC` = "CG>TC", `CG>TT` = "CG>TT", `CT>AA` = "CT>AA", `CT>AC` = "CT>AC",
    `CT>AG` = "CT>AG", `CT>GA` = "CT>GA", `CT>GC` = "CT>GC", `CT>GG` = "CT>GG",
    `CT>TA` = "CT>TA", `CT>TC` = "CT>TC", `CT>TG` = "CT>TG", `GC>AA` = "GC>AA",
    `GC>AG` = "GC>AG", `GC>AT` = "GC>AT", `GC>CA` = "GC>CA", `GC>CG` = "GC>CG",
    `GC>TA` = "GC>TA", `TA>AT` = "TA>AT", `TA>CG` = "TA>CG", `TA>CT` = "TA>CT",
    `TA>GC` = "TA>GC", `TA>GG` = "TA>GG", `TA>GT` = "TA>GT", `TC>AA` = "TC>AA",
    `TC>AG` = "TC>AG", `TC>AT` = "TC>AT", `TC>CA` = "TC>CA", `TC>CG` = "TC>CG",
    `TC>CT` = "TC>CT", `TC>GA` = "TC>GA", `TC>GG` = "TC>GG", `TC>GT` = "TC>GT",
    `TG>AA` = "TG>AA", `TG>AC` = "TG>AC", `TG>AT` = "TG>AT", `TG>CA` = "TG>CA",
    `TG>CC` = "TG>CC", `TG>CT` = "TG>CT", `TG>GA` = "TG>GA", `TG>GC` = "TG>GC",
    `TG>GT` = "TG>GT", `TT>AA` = "TT>AA", `TT>AC` = "TT>AC", `TT>AG` = "TT>AG",
    `TT>CA` = "TT>CA", `TT>CC` = "TT>CC", `TT>CG` = "TT>CG", `TT>GA` = "TT>GA",
    `TT>GC` = "TT>GC", `TT>GG` = "TT>GG")
}

cosmic_id_channel_to_type <- function(){
  c(`1:Del:C:0` = "1:Del:C", `1:Del:C:1` = "1:Del:C", `1:Del:C:2` = "1:Del:C",
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
    `5:Del:M:4` = "5:Del:M", `5:Del:M:5` = "5:Del:M")
}
