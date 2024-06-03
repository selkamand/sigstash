

fetch_as_text <- function(html, x, remove_whitespace = TRUE, remove_trailing_period = TRUE){
  txt <- html |>
    rvest::html_element(x) |>
    rvest::html_text()

  if(remove_whitespace){

    # Remove leading whitespace/newline
    txt <- sub(x=txt, pattern = "^[\n ]+", replacement = "")

    # Remove trailing whitespace/newline
    txt <- sub(x=txt, pattern = "[\n ]+$", replacement = "")
  }

  if(remove_trailing_period){
    txt <- sub(x=txt, pattern = "\\.$", replacement = "")
  }

  return(txt)
}


fetch_rna_annotations <- function(html, return_df = FALSE){
  anno <- list(
    aetiology_long = fetch_as_text(html, '#proposed-aetiology p:nth-child(2)'),
    aetiology = fetch_as_text(html, 'tbody:nth-child(5) td:nth-child(1)'),
    aetiology_support = fetch_as_text(html, 'tbody:nth-child(5) td:nth-child(2)'),
    comment = fetch_as_text(html, '#proposed-aetiology p:nth-child(4)') ,
    source_page = source_page,
    signature_version = fetch_as_text(html, '.mspage-head'),
    signature_long = fetch_as_text(html, '.mspage-subhead'),
    tissue_distribution = fetch_as_text(html, "#tissue-distribution p"),
    identification_study = fetch_as_text(html, "#checklist a"),
    first_included_in_cosmic = fetch_as_text(html, "tbody:nth-child(2) td:nth-child(2)"),
    identification_ngs_technique = fetch_as_text(html, "tbody:nth-child(3) td:nth-child(1)"),
    replicated_in_additional_studies = fetch_as_text(html, "tbody:nth-child(4) td:nth-child(2)"),
    experimental_study = fetch_as_text(html, "tbody:nth-child(6) td:nth-child(1)"),
    experimental_study_species = fetch_as_text(html, "tbody:nth-child(6) td:nth-child(2)")
  )
  anno$signature = anno$signature_long |> sub(x=_, " .*", "") |> sub(x=_, "-", ".")

  if(return_df) {
    anno <- as.data.frame(anno)
  }

  return(anno)
}

surround_by_slash <- function(x) { paste0("/", x, "/") }

source_page <- "https://cancer.sanger.ac.uk/signatures/rna-sbs/"
#source_page <- "https://cancer.sanger.ac.uk/signatures/sbs/"
#source_page <- "https://cancer.sanger.ac.uk/signatures/dbs/"
# source_page <- "https://cancer.sanger.ac.uk/signatures/id/"
# source_page <- "https://cancer.sanger.ac.uk/signatures/cn/"
# source_page <- "https://cancer.sanger.ac.uk/signatures/sv/"

link_suffix <- strsplit(source_page, split = "/", ) |>
  unlist() |>
  tail(n=2) |>
  paste0(collapse = "/") |>
  surround_by_slash()

html <- rvest::read_html(source_page)

rel_links_to_sig_pages <- html |>
  rvest::html_elements("a") |>
  rvest::html_attr("href") |>
  Filter(x = _, \(link){grepl(x=link, link_suffix, fixed = TRUE)})

urls_sig_pages <- paste0("https://cancer.sanger.ac.uk", rel_links_to_sig_pages)
urls_sig_pages


#test_html <- rvest::read_html("https://cancer.sanger.ac.uk/signatures/rna-sbs/rna-sbs1/")

# x= fetch_rna_annotations(test_html)
# x

l <- lapply(urls_sig_pages, \(url){
    sig_html <- rvest::read_html(url)
    fetch_rna_annotations(sig_html, return_df = TRUE)
  })

df_ann <- do.call("rbind", l)

df_ann[["class"]] <- "?"
df_ann[["subclass"]] <- "?"


df_ann <- df_ann |>
  dplyr::relocate(signature, class, subclass, aetiology, aetiology_long, comment,
                  signature_version, .before = dplyr::everything())

df_ann
# Write Out The data
write.csv(df_ann, file = "inst/reference_signatures/annotations/COSMIC_v3.4_RNA-SBS_annotations.csv",row.names = FALSE)

devtools::load_all(); sig_load_annotations("COSMIC_v3.4_RNA-SBS_GRCh37") |> View()


df_available <- sig_available()
datasets <- df_available[["dataset"]]

expect_gt(length(datasets), 0)

# All collections loadable without errors
for (dataset in datasets) {
  message("Testing ", dataset)
  #expect_error(sig_load_annotations(dataset), NA)

  sig <- sig_load(dataset)
  sigs <- sig[["signature"]]
  sigshared::assert_signature_annotations(sig_load_annotations(dataset), required_signatures = sigs)
}

