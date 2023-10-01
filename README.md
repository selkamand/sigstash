
<!-- README.md is generated from README.Rmd. Please edit that file -->

# sigstash

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![CRAN
status](https://www.r-pkg.org/badges/version/sigstash)](https://CRAN.R-project.org/package=sigstash)
[![R-CMD-check](https://github.com/selkamand/sigstash/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/selkamand/sigstash/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

**sigstash** makes it easy to load common mutational signature
collections into R.

To visualise or simulate decompositions using these signatures, check
out the [**sigverse**](https://github.com/selkamand/sigverse)

## Installation

You can install the development version of sigstash from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("selkamand/sigstash")
```

## Quick Start

There are only 3 commands you need to know.

1.  `sig_available` will list available signature collections

2.  `sig_load('collection_name')` will load the signature collection

3.  `sig_load_annotations('collection_name')` will load signature-level
    annotations (e.g. associated aetiologies)

### List available datasets

``` r
library(sigstash)

sig_available()
#> # A tibble: 8 × 2
#>   dataset                  description                                          
#>   <chr>                    <chr>                                                
#> 1 COSMIC_v3.3.1_SBS_GRCh38 Cancer-Related Single Base Substitution Signatures (…
#> 2 COSMIC_v3.3_DBS_GRCh38   Cancer-Related Doublet Base Substitution Signatures …
#> 3 COSMIC_v3.3_ID_GRCh37    Cancer-Related Indel Signatures (GRCh37). 83 channel…
#> 4 COSMIC_v3.3_CN_GRCh37    Cancer-Related CopyNumber Signatures (GRCh37). 48 ch…
#> 5 COSMIC_v3.3.1_SBS_GRCh37 Cancer-Related Single Base Substitution Signatures (…
#> 6 COSMIC_v3.3_DBS_GRCh37   Cancer-Related Doublet Base Substitution Signatures …
#> 7 COSMIC_v2_SBS_GRCh38     An older set of cancer-related Single Base Substitut…
#> 8 COSMIC_v2_SBS_GRCh37     An older set of cancer-related Single Base Substitut…
```

### Load and Visualise datasets

``` r
cosmic_signatures = sig_load('COSMIC_v3.3.1_SBS_GRCh38')
```

### Load signature metadata

Some signature collections will include signature level annotations for
example describing aetiologies typically associated with signatures.
These can be loaded using `sig_load_annotations`

``` r
sig_load_annotations('COSMIC_v3.3.1_SBS_GRCh38')
#> # A tibble: 79 × 14
#>    signa…¹ class subcl…² aetio…³ aetio…⁴ comment Signa…⁵ Ident…⁶ Ident…⁷ Exper…⁸
#>    <chr>   <chr> <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>  
#>  1 SBS1    cloc… clock-… " Spon… An end… "Signa… Mutati… Nik-Za… https:… ""     
#>  2 SBS2    cyti… cytidi… " APOB… Attrib… "SBS2 … Mutati… Nik-Za… https:… "Chan …
#>  3 SBS3    dysf… HR      " HR d… Defect… "SBS3 … Mutati… Nik-Za… https:… "Zámbo…
#>  4 SBS4    envi… tobacco " Agin… Associ… "Altho… Mutati… Alexan… https:… "Nik-Z…
#>  5 SBS5    cloc… clock-… " Agin… Unknow… "SBS5 … Mutati… Alexan… https:… ""     
#>  6 SBS6    dysf… MMR     " MMR … SBS6 i… "SBS6 … Mutati… Alexan… https:… "Meier…
#>  7 SBS7a   envi… UV      " UV l… SBS7a/… ""      Mutati… Haywar… https:… "Nik-Z…
#>  8 SBS7b   envi… UV      " UV l… SBS7a/… ""      Mutati… Haywar… https:… "Nik-Z…
#>  9 SBS7c   envi… UV      " UV l… SBS7a/… ""      Mutati… Saini … https:… ""     
#> 10 SBS7d   envi… UV      " UV l… SBS7a/… ""      Mutati… Saini … https:… ""     
#> # … with 69 more rows, 4 more variables: ExperimentalValidationURL <chr>,
#> #   ProposedAetiologySupport <chr>, ValidatedInOrthogonalTechniques <chr>,
#> #   source.page <chr>, and abbreviated variable names ¹​signature, ²​subclass,
#> #   ³​aetiology, ⁴​aetiology_long, ⁵​SignatureVersion, ⁶​IdentificationStudy,
#> #   ⁷​IdentificationStudyURL, ⁸​ExperimentalValidationStudy
```
