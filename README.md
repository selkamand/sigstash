
<!-- README.md is generated from README.Rmd. Please edit that file -->

# sigstash

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![CRAN
status](https://www.r-pkg.org/badges/version/sigstash)](https://CRAN.R-project.org/package=sigstash)
[![R-CMD-check](https://github.com/selkamand/sigstash/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/selkamand/sigstash/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/selkamand/sigstash/branch/main/graph/badge.svg)](https://app.codecov.io/gh/selkamand/sigstash?branch=main)
![GitHub Issues or Pull
Requests](https://img.shields.io/github/issues-closed/selkamand/sigstash)
[![code
size](https://img.shields.io/github/languages/code-size/selkamand/sigstash.svg)](https://github.com/selkamand/sigstash)
![GitHub last
commit](https://img.shields.io/github/last-commit/selkamand/sigstash)
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
#> # A tibble: 25 × 2
#>    dataset                    description                                       
#>    <chr>                      <chr>                                             
#>  1 REFSIG_v2.03_SBS           Cancer-Related Single Base Substitution Signature…
#>  2 REFSIG_v2.03_DBS           Cancer-Related Doublet Base Substitution Signatur…
#>  3 COSMIC_v3.4_SBS_GRCh38     Cancer-Related Single Base Substitution Signature…
#>  4 COSMIC_v3.4_DBS_GRCh38     Cancer-Related Doublet Base Substitution Signatur…
#>  5 COSMIC_v3.4_SV_GRCh38      Cancer-Related Structural Variant  Signatures (GR…
#>  6 COSMIC_v3.4_CN_GRCh37      Cancer-Related CopyNumber Signatures (GRCh38). 48…
#>  7 COSMIC_v3.4_DBS_GRCh37     Cancer-Related Doublet Base Substitution Signatur…
#>  8 COSMIC_v3.4_ID_GRCh37      Cancer-Related Indel Signatures (GRCh37). 83 chan…
#>  9 COSMIC_v3.4_SBS_GRCh37     Cancer-Related Single Base Substitution Signature…
#> 10 COSMIC_v3.4_RNA-SBS_GRCh37 Cancer-Related Single Base Substitutions in RNA s…
#> # ℹ 15 more rows
```

### Load and Visualise datasets

``` r
# Load Collection
cosmic_signatures <- sig_load("COSMIC_v3.3.1_SBS_GRCh38")

# Access a specific signature
cosmic_signatures[["SBS1"]]
#> # A tibble: 96 × 3
#>    type  channel fraction
#>    <chr> <chr>      <dbl>
#>  1 C>A   A[C>A]A 0.000876
#>  2 C>A   A[C>A]C 0.00222 
#>  3 C>A   A[C>A]G 0.000180
#>  4 C>A   A[C>A]T 0.00127 
#>  5 C>G   A[C>G]A 0.00184 
#>  6 C>G   A[C>G]C 0.00119 
#>  7 C>G   A[C>G]G 0.000117
#>  8 C>G   A[C>G]T 0.00113 
#>  9 C>T   A[C>T]A 0.0247  
#> 10 C>T   A[C>T]C 0.00615 
#> # ℹ 86 more rows
```

### Load signature metadata

Some signature collections will include signature level annotations for
example describing aetiologies typically associated with signatures.
These can be loaded using `sig_load_annotations`

``` r
sig_load_annotations("COSMIC_v3.3.1_SBS_GRCh38")
#> # A tibble: 79 × 14
#>    signature class    subclass aetiology aetiology_long comment SignatureVersion
#>    <chr>     <chr>    <chr>    <chr>     <chr>          <chr>   <chr>           
#>  1 SBS1      clock-l… clock-l… " Sponta… An endogenous… "Signa… Mutational Sign…
#>  2 SBS2      cytidin… cytidin… " APOBEC… Attributed to… "SBS2 … Mutational Sign…
#>  3 SBS3      dysfunc… HR       " HR def… Defective hom… "SBS3 … Mutational Sign…
#>  4 SBS4      environ… tobacco  " Aging … Associated wi… "Altho… Mutational Sign…
#>  5 SBS5      clock-l… clock-l… " Aging … Unknown. SBS5… "SBS5 … Mutational Sign…
#>  6 SBS6      dysfunc… MMR      " MMR de… SBS6 is assoc… "SBS6 … Mutational Sign…
#>  7 SBS7a     environ… UV       " UV lig… SBS7a/SBS7b/S… ""      Mutational Sign…
#>  8 SBS7b     environ… UV       " UV lig… SBS7a/SBS7b/S… ""      Mutational Sign…
#>  9 SBS7c     environ… UV       " UV lig… SBS7a/SBS7b/S… ""      Mutational Sign…
#> 10 SBS7d     environ… UV       " UV lig… SBS7a/SBS7b/S… ""      Mutational Sign…
#> # ℹ 69 more rows
#> # ℹ 7 more variables: IdentificationStudy <chr>, IdentificationStudyURL <chr>,
#> #   ExperimentalValidationStudy <chr>, ExperimentalValidationURL <chr>,
#> #   ProposedAetiologySupport <chr>, ValidatedInOrthogonalTechniques <chr>,
#> #   source.page <chr>
```
