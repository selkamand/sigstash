devtools::load_all()
path_inst <- system.file(package = "sigstash")
path_spreadsheet <-  paste0(path_inst, "/original_files/signal/SupplementaryTables.xlsx")

# Originally downloaded from https://github.com/Nik-Zainal-Group/signature.tools.lib/tree/master/data/RefSigSBS_v2.03
#df_refsigs <- utils::read.csv(paste0(path_inst, "/original_files/signal/RefSig_SBS_v2.03.tsv"), sep = "\t")

# Originally downloaded from https://www.science.org/doi/10.1126/science.abl9283
df_refsigs_sbs <- readxl::read_excel(path_spreadsheet,sheet = "Table S21", )



df_refsigs_sbs_long <- df_refsigs_sbs |>
  dplyr::rename(channel = mutationClass) |>
  tidyr::pivot_longer(cols = -channel, names_to = "signature", values_to = "fraction", values_transform = as.numeric)

df_refsigs_sbs_long[["type"]] <- sigstash::sig_convert_channel2type(df_refsigs_sbs_long[["channel"]],sigclass = "SBS96")
df_refsigs_sbs_long

df_refsigs_long_sorted <- df_refsigs_sbs_long |>
  dplyr::relocate(signature, type, channel, fraction) |>
  dplyr::arrange(readr::parse_number(signature))

df_refsigs_long_sorted

# df_refsigs_long_sorted |>
#   write.csv("inst/reference_signatures/REFSIG_v2.03_SBS.csv")


# Annotations -------------------------------------------------------------


# Annotations
df_refsigs_sbs_annotations <- readxl::read_excel(path_spreadsheet,sheet = "Table S19")
df_refsigs_sbs_annotations <- df_refsigs_sbs_annotations |>
  dplyr::rename(
    signature_fullname = signature,
    aetiology = Aetiology,
    )

df_refsigs_sbs_annotations <- df_refsigs_sbs_annotations |>
  dplyr::mutate(
    signature = stringr::str_replace(signature_fullname, ".* ", "")
  )

df_refsigs_sbs_annotations <- df_refsigs_sbs_annotations |>
  dplyr::select(signature_fullname, signature, QC, aetiology, comment = notes, refsigBreakDown, TSB, RSB, COSMICv3.2, RefSigv1) |>
  dplyr::mutate(subclass = dplyr::case_when(
      aetiology == "Smoking" ~ "tobacco",
      aetiology == "HR deficiency" ~ "HR",
      aetiology == "MMR deficiency" ~ "MMR",
      aetiology == "APOBEC" ~ "cytidine_deaminases",
      aetiology == "AID" ~ "cytidine_deaminases",
      aetiology == "Deamination" ~ "cytosine_deamination",
      aetiology == "Platinum-based theraphy" ~ "chemotherapy_platinum",
      aetiology == "Thiopurine" ~ "chemotherapy_thiopurine",
      aetiology == "Azathioprine treatment" ~ "immunosuppression",
      aetiology == "Duocarmycin" ~ "chemotherapy_nitrogen_mustards",
      aetiology == "temozolomide/1,2-DMH" ~ "triazenes",
      aetiology == "Aflatoxin" ~ "aflatoxin",
      aetiology == "Aristolochic acid" ~ "aristolochic_acid",
      aetiology == "POLE dysregulation" ~ "polymerase_mutations",
      aetiology == "MMR+POLE deficiency" ~ "MMR",
      aetiology == "MMR+POLD deficiency" ~ "MMR",
      aetiology == "MMR deficiency (PMS2)" ~ "MMR",
      aetiology == "NTHL1 deficiency" ~ "BER",
      aetiology == "ROS damage/MUTYH" ~ "BER",
      aetiology == "UV light" ~ "UV",
      aetiology == "Age of diagnosis" ~ "clock-like",
      aetiology == "Colibactin" ~ "colibactin",
      aetiology == "Possible sequencing artefact" ~ "sequencing_artefact",
      aetiology == "Lymphoma" ~ "polymerase_mutations",
      aetiology == "Treatment" ~ "unknown",
      aetiology == "Only found in ICGC-Liver" ~ "unknown",
      aetiology == "Unknown" ~ "unknown",
      .default = NA
    ))

df_aetiology_classes <- sigshared::sig_aetiology_classes()

df_refsigs_sbs_annotations_withclass <- df_refsigs_sbs_annotations |>
  dplyr::left_join(df_aetiology_classes, by = "subclass")

df_refsigs_sbs_annotations_withclass |>
  dplyr::filter(is.na(class))

df_refsigs_sbs_annotations_withclass |>
  dplyr::mutate(aetiology_long = aetiology) |>
  dplyr::relocate(signature, QC, class, subclass, aetiology, aetiology_long, comment) |>
  utils::write.csv("inst/reference_signatures/annotations/REFSIG_v2.03_SBS_annotations.csv", row.names = FALSE)

