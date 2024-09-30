#' Write a Signature Collection to a File
#'
#' This function writes a signature collection, in sigstash format, to a specified file. The collection can be saved either as an RDS object, a tidy CSV, or a wide CSV.
#'
#' @param signatures A signature collection in sigstash format. This should be a list of signatures.
#' @param format The format in which the signatures should be saved. Options are:
#' \describe{
#'   \item{\strong{"rds"}}{Serializes the list in RDS format.}
#'   \item{\strong{"csv_tidy"}}{Converts the collection to a tidy (long) data frame and stores it as a CSV file.}
#'   \item{\strong{"csv_wide"}}{Converts the collection to a wide format (with channels as row names) and stores it as a CSV file.}
#' }
#' @param filepath The file path where the signature collection should be saved.
#' @return Invisibly returns \code{NULL}. This function is called for its side effects of writing data to a file.
#' @export
#'
#' @examples
#' \dontrun{
#' signature_collection <- sigstash::sig_load("COSMIC_v3.4_SBS_GRCh38")
#' sig_write_collection(signature_collection, filepath = "signatures.rds", format = "rds")
#' }
sig_write_signatures <- function(signatures, filepath, format = c("rds", "csv_tidy", "csv_wide")) {
  sigshared::assert_signature_collection(signatures)
  assertions::assert_string(filepath)
  format <- rlang::arg_match(format)

  if(format == "rds") {
    saveRDS(signatures, file = filepath, compress = TRUE)
  } else if(format == "csv_tidy") {
    df_tidy <- sig_collection_reformat_list_to_tidy(signatures)
    utils::write.csv(df_tidy, file = filepath, row.names = FALSE)
  } else if(format == "csv_wide") {
    df_wide <- sig_collection_reformat_list_to_wide(signatures)
    utils::write.csv(df_wide, file = filepath, row.names = TRUE)
  } else {
    stop("Serialization of signature collection to format: ", format, " has not yet been implemented")
  }

  return(invisible(NULL))
}

#' Read a Signature Collection from a File
#'
#' This function reads a signature collection from a specified file, which can be either an RDS file or a tidy CSV. The collection is returned as a list of signatures in sigstash format.
#'
#' @param filepath The file path from which the signature collection should be read.
#' @param format The format of the file to be read. Options are:
#' \describe{
#'   \item{\strong{"rds"}}{Reads the collection from an RDS file.}
#'   \item{\strong{"csv_tidy"}}{Reads the collection from a tidy (long) CSV file and converts it to a list of signatures.}
#' }
#' @param collection_name A string specifying the signature collection name to be added as the 'collection_name' attribute
#' of the parsed signature object.
#' If NULL, name is guessed from the filename (before the first '.'). Ignored if `format = 'rds'`.
#' @return A list of signatures in sigstash format.
#' @export
#'
#' @note Wide-format CSVs cannot be read back into sigverse format because they lack the necessary channel type information for accurate reconstruction (channel type is not included).
#'
#' @examples
#' \dontrun{
#' signatures <- sig_read_collection(filepath = "signatures.rds", format = "rds")
#' }
sig_read_signatures <- function(filepath, format = c("rds", "csv_tidy"), collection_name=NULL) {

  # Assertions
  assertions::assert_file_exists(filepath)
  if(!is.null(collection_name)) assertions::assert_string(collection_name)

  format <- rlang::arg_match(format)

  # If RDS serialisation
  if(format == "rds") {
    signatures <- readRDS(filepath)
    return(signatures)
  }
  else if(format == "csv_tidy") {
    df_tidy <- utils::read.csv(file = filepath, header = TRUE)
    signatures <- sig_collection_reformat_tidy_to_list(df_tidy)

    # Add attributes to signature collection
    if(is.null(collection_name)) {
      collection_name <- sub(x=basename(filepath), pattern = "\\..*", "")
      }
    signatures <- add_collection_attributes(signatures, name = collection_name, format = "sigstash")

    return(signatures)
  }
  else {
    stop("Reading signature collections from format: ", format, " has not been implemented yet")
  }
}
