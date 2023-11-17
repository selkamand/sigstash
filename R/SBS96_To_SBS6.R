
#' Reduce Signature Complexity
#'
#' Reduce complex signatures into simpler forms by converting 'channel' level signatures into 'type' level signatures Takes complex signatures like SBS96 and turns them into SBS6 signatures.
#'
#' @param signatures a sigverse signature collection
#'
#' @return a list conforming to sigverse signature collection format
#'
#' @export
sig_reduce_complexity_using_type <- function(signatures){
   lapply(signatures, \(sig){
     df = stats::aggregate(fraction ~ type, data = sig, sum)
     df[['channel']] <- df[['type']]
     df <- tibble::tibble(df[c('type', 'channel', 'fraction')])
     return(df)
   })
}

