devtools::load_all()


experimental_datasets <- sig_available() |> dplyr::filter(grepl(x=dataset, "EXPERIMENTAL")) |> dplyr::pull(dataset)
experimental_datasets

lapply(experimental_datasets, \(dataset){
  df <- dplyr::tibble(
    signature = sig_load(dataset) |> names(),
    aetiology =  sub(x=signature,pattern = "_.*?$", ""),
    class = "?",
    subclass = "?"
  )

  annotation_filepath <- dataset |>
    strsplit("_") |>
    lapply(\(x) { head(x, n=5)  |> paste0(collapse = "_") })  |>
    unlist() |>
    paste0("inst/reference_signatures/", .=_, "_annotations.csv")

  readr::write_csv(df, annotation_filepath)
})



experimental_datasets |>
  strsplit("_") |>
  lapply(\(x) { head(x, n=5)  |> paste0(collapse = "_") })  |>
  unlist() |>
  paste0("inst/reference_signatures/", .=_, "_annotations.csv") |>
  clipr::write_clip()

df

