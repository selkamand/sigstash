# Load datasets.csv
load_datasets <- function() {
  path_datasets <- system.file("datasets.csv", package = "sigstash")
  df_datasets <- utils::read.csv(path_datasets, header = TRUE)
  return(df_datasets)
}


#' List available Signature Collection
#'
#' Lists all available signature collections available in sigstash
#'
#' @return a data.frame describing available datasets
#' @export
#'
#' @inherit sig_load examples
sig_available <- function() {
  df_datasets <- load_datasets()
  df_datasets <- df_datasets[, c("dataset", "description")]

  df_datasets <- tibble::tibble(df_datasets)
  return(df_datasets)
}


#' Load a Signature Collection
#'
#' @param dataset a valid signature collection to load. Run [sig_available()] to list available datasets.
#' @param format what format should we return signature collection in (string)
#'
#' @return The specified signature collection, either as a list conforming to sigverse signature collection format, or, if `return_df = TRUE` as a data.frame.
#' @export
#'
#' @examples
#'
#' # List available datasets
#' sig_available()
#'
#' # Load available datasets
#' sig_load("COSMIC_v3.3.1_SBS_GRCh38")
#'
#' @details
#' Use the format \strong{format} argument to change return format:
#'
#' |                   |                                                                                               |
#' |-------------------|-----------------------------------------------------------------------------------------------|
#' | \strong{sigstash} | Signatures returned as a list of dataframes where each data.frame is a signature              |
#' | \strong{tidy}     | Signatures returned as a tidy 4-column dataframe with signature, type, channel, fraction      |
#' | \strong{sigminer} | Signatures returned as a single dataframe where columns are samples and rows are channels. Compatible with sigminer |
sig_load <- function(dataset, format = c("sigstash", "tidy", "sigminer")) {
  # Assertions
  assertions::assert_string(dataset)
  format <- rlang::arg_match(format)

  df_datasets <- load_datasets()
  valid_datasets <- df_datasets[["dataset"]]

  assertions::assert_subset(dataset, valid_datasets, msg = "`{arg_value}` is not a valid sigstash dataset. See {.code sig_available()} for a list of datasets")

  path <- df_datasets[["path_signature"]][match(dataset, df_datasets[["dataset"]])]
  path <- system.file(path, package = "sigstash")

  df_data <- utils::read.csv(path, header = TRUE)

  ls_data <- sig_collection_reformat_tidy_to_list(df_data)

  if (format == "sigstash") {
    return(ls_data)
  }

  if (format == "tidy") {
    df_data <- sig_collection_reformat_list_to_tidy(ls_data)
    return(df_data)
  }

  if (format == "sigminer") {
    df_data <- sig_collection_to_sigminer(ls_data)
    return(df_data)
  }

  stop("Unexpected value of `format` argument: ", paste0(format, collapse = ","))

  return(invisible(NULL))
}


#' Load Signature Annotations
#'
#' Often signature collections contain annotations describing, for example, aetiology of each signature in collections
#'
#' @param dataset a valid signature collection to load. Run [sig_available()] to list available datasets.
#'
#' @return A data.frame containing all signature-level annotations of a given signature collection.
#' @export
#'
#' @examples
#'
#' # List available datasets
#' sig_available()
#'
#' # Load available datasets
#' signature_collection <- sig_load("COSMIC_v3.3.1_SBS_GRCh38")
#'
#' # Load associated annotations
#' signature_annotations <- sig_load_annotations("COSMIC_v3.3.1_SBS_GRCh38")
#'
sig_load_annotations <- function(dataset) {
  # Assertions
  assertions::assert_string(dataset)


  df_datasets <- load_datasets()
  valid_datasets <- df_datasets[["dataset"]]

  assertions::assert_subset(dataset, valid_datasets, msg = "`{arg_value}` is not a valid sigstash dataset. See {.code sig_available()} for a list of datasets")

  # Read in signature collection
  annotation_file_no_extension <- df_datasets[["annotations"]][match(dataset, df_datasets[["dataset"]])]
  path <- system.file("reference_signatures/annotations/", paste0(annotation_file_no_extension, ".csv"), package = "sigstash")

  assertions::assert_file_exists(path, msg = "No signature annotation file found for dataset: {.strong {dataset}}")
  df_data <- utils::read.csv(path, header = TRUE)

  # Return annotation data.frame
  df_data <- tibble::tibble(df_data)
  return(df_data)
}
