#' Signature Collection Set Operations
#'
#' Create custom signature collections using [sig_subset()] and [sig_merge()].
#' Note sig_merge concatenates multiple signature collections together,
#' and should not be confused with [sigstats::sig_combine()]
#' which generates signature models by mathematically combining signatures from a single collection.
#'
#' @param signatures A sigverse signature collection. A named list of sigstash signature data.frames. See \href{https://github.com/selkamand/sigshared?tab=readme-ov-file#sigverse-data-types}{sigshared readme} for details.
#' @param names A character vector describing the signatures to subset.
#'
#' @return A sigverse signature collection (named list of signature data.frames) containing only the signatures specified by the \code{names} argument
#' @export
#'
#' @examples
#' # Load a signature collection
#' signature_collection <- sigstash::sig_load("COSMIC_v3.4_SBS_GRCh38")
#'
#' # Subset collection to include only 2 signatures
#' custom_signature_collection <- sig_subset(signature_collection, c("SBS3", "SBS8"))
#'
#' # Load a second signature collection
#' signature_collection2 <- sigstash::sig_load("COSMIC_EXPERIMENTAL_V1_SBS_HUMAN_UNFILTERED")
#'
#' # Subset collection to include only 1 signature
#' custom_signature_collection2 <- sig_subset(signature_collection2, "temozolomide_8377cb5da514")
#'
#' # Merge our 2 collections together
#' sig_merge(custom_signature_collection, custom_signature_collection2)
sig_subset <- function(signatures, names) {
  # Assertions
  sigshared::assert_signature_collection(signatures)
  assertions::assert_names_include(signatures, names)

  # Return Subset
  signatures[names]
}


#' Merge Signature Collections
#'
#' @param ... Signature collections to merge.
#'
#' @return A merged list of signature collections. See \href{https://github.com/selkamand/sigshared?tab=readme-ov-file#sigverse-data-types}{sigshared readme} for details.
#' @export
#'
#' @inherit sig_subset examples description
#'
sig_merge <- function(...) {
  # Combine all lists into one
  all_signatures <- list(...)

  # Flatten the list of lists into a single list
  combined_signatures <- do.call(c, all_signatures)

  # Check for duplicate names
  signature_names <- names(combined_signatures)
  duplicates <- signature_names[duplicated(signature_names)]
  assertions::assert_no_duplicates(signature_names, msg = "Duplicate signatures found [{duplicates}]. Please rename duplicate signatures before merging")

  # Return the combined list
  combined_signatures
}
