devtools::load_all()

sig_available()

dataset_old ="COSMIC_v3.3_CN_GRCh37"
dataset_new = "COSMIC_v3.4_CN_GRCh37"
ann1=sig_load_annotations(dataset_old)
ann2=sig_load_annotations(dataset_new)
ann2 |>
  dplyr::count(class, sort = TRUE, subclass)

port_cols <- ann1 |>
  dplyr::select(signature, class, subclass) |>
  dplyr::rename(class = class, subclass = subclass)



ann2 |>
  dplyr::select(-class, -subclass) |>
  dplyr::left_join(port_cols, by = "signature") |>
  dplyr::relocate(class, subclass, .after = signature)  |>
  dplyr::select(class, subclass) |>
  clipr::write_clip()





sig_load_annotations("COSMIC_v3.4_SBS_GRCh38") |>
  dplyr::count(class, subclass, sort = TRUE) |> View()
