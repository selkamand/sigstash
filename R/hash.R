#' Check if a signature collection matches a known dataset
#'
#' The `sig_identify_collection` function checks whether a given signature collection
#' matches any precomputed dataset by computing its MD5 checksum and comparing it
#' with a reference set of checksums. The reference checksums are loaded from
#' precomputed data. The function can return either the dataset name or the
#' dataset name along with the format (e.g., 'tidy', 'sigstash', 'sigminer').
#'
#' @param signatures A signature dataset (in sigstash, tidy, or sigminer format) to be checked.
#'   The signatures will be sorted before the MD5 checksum is computed.
#' @param return A character string specifying the type of output. Should be one of:
#'   - `"name"`: Returns just the dataset name.
#'   - `"name_plus_format"`: Returns both the dataset name and format in a single string.
#'   Default is `"name"`.
#'
#' @return
#' If the signature collection matches a known dataset, the function returns either
#' the dataset name or the dataset name with the format. If no match is found,
#' it returns `"Uncertain"`.
#'
#' @examples
#' sig_identify_collection(sig_load("COSMIC_v3.4_SBS_GRCh38"), return = "name")
#'
#' sig_identify_collection(sig_load("COSMIC_v3.4_SBS_GRCh38"), return = "name_plus_format")
#'
#' @export
sig_identify_collection <- function(signatures, return = c("name", "name_plus_format", "md5")) {

  # Ensure 'return' is one of the allowed options
  return <- rlang::arg_match(return)

  # Load the precomputed MD5 sums for all known datasets
  md5_to_collection_map <- load_precomputed_md5s()

  # Compute the MD5 checksum of the input signature collection
  md5sum <- get_md5sum(signatures)

  if(return == "md5") {return(md5sum)}

  # If the MD5 checksum is not found in the precomputed collection, return "Uncertain"
  if (!md5sum %in% names(md5_to_collection_map)) {
    return("Uncertain")
  }

  # Retrieve the description corresponding to the MD5 sum
  description <- md5_to_collection_map[[md5sum]]

  # Return the appropriate value based on the 'return' argument
  if (return == "name") {
    return(description[["dataset"]])
  } else if (return == "name_plus_format") {
    return(description[["unique_name"]])
  } else {
    stop("Unexpected value of return: ", as.character(return))
  }
}

precompute_and_save_md5s <- function(){
  requireNamespace("here", quietly = TRUE)
  filepath = here::here("inst/reference_signatures/md5sums.Rds")
  message("Writing md5s to [", filepath, "]")
  l <- compute_md5_for_all_collections()
  saveRDS(l, file = filepath)
}


#' Compute MD5 sums for all datasets and formats
#'
#' Computes MD5 checksums for each dataset in available signature collections,
#' across formats ('tidy', 'sigstash', 'sigminer'). Before calculating the
#' checksums, the signatures are sorted internally, so their prior sorting
#' does not affect the result.
#' Numeric values are also rounded to 5 decimal places to avoid any OS-specific
#' issues.
#'
#' @return
#' A named list where each MD5 sum is a key, and the value is a list containing:
#' - `unique_name`: A string describing the dataset and format.
#' - `dataset`: The dataset name.
#' - `format`: The format of the dataset ('tidy', 'sigstash', or 'sigminer').
#'
#' @note This function is internal and will not be exported.
#'
#' @keywords internal
compute_md5_for_all_collections <- function(){
  available <- sig_available()
  datasets <- available[['dataset']]

  # Initialize an empty list to store all keys
  all_keys <- list()

  # Loop over each dataset
  for (dataset_name in datasets) {
    # Compute MD5 sums for each format
    md5sum_sigstash <- get_md5sum(sig_load(dataset_name, format = "sigstash"))
    md5sum_sigminer <- get_md5sum(sig_load(dataset_name, format = "sigminer"))
    md5sum_tidy     <- get_md5sum(sig_load(dataset_name, format = "tidy"))

    # Create entries for each MD5 sum
    all_keys[[md5sum_tidy]] <- list(
      unique_name = paste(dataset_name, "(format: tidy)"),
      dataset = dataset_name,
      format = "tidy"
    )

    all_keys[[md5sum_sigstash]] <- list(
      unique_name = paste(dataset_name, "(format: sigstash)"),
      dataset = dataset_name,
      format = "sigstash"
    )

    all_keys[[md5sum_sigminer]] <- list(
      unique_name = paste(dataset_name, "(format: sigminer)"),
      dataset = dataset_name,
      format = "sigminer"
    )
  }

  # Return the lookup list
  return(all_keys)
}



load_precomputed_md5s <- function(){
  path = system.file("reference_signatures/md5sums.Rds", package = "sigstash")
  assertions::assert_file_exists(path)
  readRDS(path)
}

load_precomputed_md5s_as_df <- function(){
  l <- load_precomputed_md5s()
  tibble::tibble(
    md5 = names(l),
    unique_name = vapply(l, \(l2){l2$unique_name}, FUN.VALUE = character(1))
    )
}


round_fractions <- function(dataset, digits = 5){
  dataset[] <- lapply(dataset, function(x) {
    if(is.numeric(x)) round(x, digits) else x
  })
}

get_md5sum <- function(dataset){
  if(is.data.frame(dataset)){
    dataset <- sort_signatures_dataframe(dataset)
    dataset <- round_fractions(dataset)
  }
  else if(is.list(dataset)){
    dataset <- sort_signatures_list(dataset)
    dataset <- lapply(dataset, round_fractions)
  }
  else
    stop("Unexpected dataset class")

  digest::digest(dataset, algo = "md5")
}
sort_signatures_list <- function(l){
  l <- l[li_order(names(l))]
  lapply(l, FUN = \(df) { sort_signatures_dataframe(df)})
}

sort_signatures_dataframe <- function(df){
  if ("channel" %in% colnames(df)) {
    return(df[li_order(df[["channel"]]),])
  }
  else{
    return(df[li_order(row.names(df)), ])
  }
}

# Locale independent order function
li_order <- function(expr){
  # Save the current locale setting
  original_locale <- Sys.getlocale("LC_COLLATE")

  # Ensure the original locale is restored even if an error occurs
  on.exit(Sys.setlocale("LC_COLLATE", original_locale), add = TRUE)

  # Set the locale to "C" for consistent ordering
  Sys.setlocale("LC_COLLATE", "C")

  # Perform Ordering
  order(expr)
}
