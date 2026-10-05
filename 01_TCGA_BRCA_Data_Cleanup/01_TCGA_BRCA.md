---
title: "TCGA-BRCA: RNA-seq + clinical data cleanup (ER+ MBC vs FBC cohort)"
author: "Junkai Wang"
date: "Original: 7/1/2020; revised: 9/4/2026"
output:
  html_document:
    keep_md: true
    toc: true
editor_options: 
  chunk_output_type: console
---


``` r
knitr::opts_chunk$set(message = FALSE, warning = FALSE)
knitr::opts_knit$set(root.dir = ".")

raw_dir <- "Raw_Data"
out_dir <- "Output"
dir.create(raw_dir, showWarnings = FALSE)
dir.create(out_dir, showWarnings = FALSE)

f_se     <- file.path(raw_dir, "TCGA_BRCA_SE_raw.rds")
f_biotab <- file.path(raw_dir, "TCGA_BRCA_biotab.rds")

library(TCGAbiolinks)
library(SummarizedExperiment)
```

```
## Loading required package: MatrixGenerics
```

```
## Loading required package: matrixStats
```

```
## 
## Attaching package: 'MatrixGenerics'
```

```
## The following objects are masked from 'package:matrixStats':
## 
##     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
##     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
##     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
##     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
##     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
##     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
##     colWeightedMeans, colWeightedMedians, colWeightedSds,
##     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
##     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
##     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
##     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
##     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
##     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
##     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
##     rowWeightedSds, rowWeightedVars
```

```
## Loading required package: GenomicRanges
```

```
## Loading required package: stats4
```

```
## Loading required package: BiocGenerics
```

```
## Loading required package: generics
```

```
## 
## Attaching package: 'generics'
```

```
## The following objects are masked from 'package:base':
## 
##     as.difftime, as.factor, as.ordered, intersect, is.element, setdiff,
##     setequal, union
```

```
## 
## Attaching package: 'BiocGenerics'
```

```
## The following objects are masked from 'package:stats':
## 
##     IQR, mad, sd, var, xtabs
```

```
## The following objects are masked from 'package:base':
## 
##     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
##     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
##     get, grep, grepl, is.unsorted, lapply, Map, mapply, match, mget,
##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
##     rbind, Reduce, rownames, sapply, saveRDS, table, tapply, unique,
##     unsplit, which.max, which.min
```

```
## Loading required package: S4Vectors
```

```
## 
## Attaching package: 'S4Vectors'
```

```
## The following object is masked from 'package:utils':
## 
##     findMatches
```

```
## The following objects are masked from 'package:base':
## 
##     expand.grid, I, unname
```

```
## Loading required package: IRanges
```

```
## Loading required package: Seqinfo
```

```
## Loading required package: Biobase
```

```
## Welcome to Bioconductor
## 
##     Vignettes contain introductory material; view with
##     'browseVignettes()'. To cite Bioconductor, see
##     'citation("Biobase")', and for packages 'citation("pkgname")'.
```

```
## 
## Attaching package: 'Biobase'
```

```
## The following object is masked from 'package:MatrixGenerics':
## 
##     rowMedians
```

```
## The following objects are masked from 'package:matrixStats':
## 
##     anyMissing, rowMedians
```

``` r
library(data.table)
```

```
## 
## Attaching package: 'data.table'
```

```
## The following object is masked from 'package:SummarizedExperiment':
## 
##     shift
```

```
## The following object is masked from 'package:GenomicRanges':
## 
##     shift
```

```
## The following object is masked from 'package:IRanges':
## 
##     shift
```

```
## The following objects are masked from 'package:S4Vectors':
## 
##     first, second
```

``` r
library(tidyverse)
```

```
## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
## ✔ dplyr     1.2.1     ✔ readr     2.2.0
## ✔ forcats   1.0.1     ✔ stringr   1.6.0
## ✔ ggplot2   4.0.3     ✔ tibble    3.3.1
## ✔ lubridate 1.9.5     ✔ tidyr     1.3.2
## ✔ purrr     1.2.2
```

```
## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
## ✖ lubridate::%within%() masks IRanges::%within%()
## ✖ dplyr::between()      masks data.table::between()
## ✖ dplyr::collapse()     masks IRanges::collapse()
## ✖ dplyr::combine()      masks Biobase::combine(), BiocGenerics::combine()
## ✖ dplyr::count()        masks matrixStats::count()
## ✖ dplyr::desc()         masks IRanges::desc()
## ✖ tidyr::expand()       masks S4Vectors::expand()
## ✖ dplyr::filter()       masks stats::filter()
## ✖ dplyr::first()        masks data.table::first(), S4Vectors::first()
## ✖ lubridate::hour()     masks data.table::hour()
## ✖ lubridate::isoweek()  masks data.table::isoweek()
## ✖ lubridate::isoyear()  masks data.table::isoyear()
## ✖ dplyr::lag()          masks stats::lag()
## ✖ dplyr::last()         masks data.table::last()
## ✖ lubridate::mday()     masks data.table::mday()
## ✖ lubridate::minute()   masks data.table::minute()
## ✖ lubridate::month()    masks data.table::month()
## ✖ ggplot2::Position()   masks BiocGenerics::Position(), base::Position()
## ✖ lubridate::quarter()  masks data.table::quarter()
## ✖ purrr::reduce()       masks GenomicRanges::reduce(), IRanges::reduce()
## ✖ dplyr::rename()       masks S4Vectors::rename()
## ✖ lubridate::second()   masks data.table::second(), S4Vectors::second()
## ✖ lubridate::second<-() masks S4Vectors::second<-()
## ✖ dplyr::slice()        masks IRanges::slice()
## ✖ purrr::transpose()    masks data.table::transpose()
## ✖ lubridate::wday()     masks data.table::wday()
## ✖ lubridate::week()     masks data.table::week()
## ✖ lubridate::yday()     masks data.table::yday()
## ✖ lubridate::year()     masks data.table::year()
## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors
```

## 1. Download RNA-seq expression (STAR counts, harmonized)


``` r
query_exp <- GDCquery(
  project       = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type     = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)
GDCdownload(query_exp, files.per.chunk = 50, directory = file.path(raw_dir, "GDCdata"))
se <- GDCprepare(query_exp, directory = file.path(raw_dir, "GDCdata"))
assays(se) <- assays(se)[c("unstranded", "tpm_unstrand")]
saveRDS(se, f_se)                                 # checkpoint
```


``` r
se <- readRDS(f_se)
dim(se)
```

```
## [1] 60660  1231
```

``` r
### ~60,660 genes x ~1,230 samples (all sample types) ###
```

## 2. Clinical data A: indexed clinical (sex, age, stage, survival)


``` r
clin <- GDCquery_clinic(project = "TCGA-BRCA", type = "clinical")
dim(clin)
```

```
## [1] 1098  100
```

``` r
clin_core <- clin %>%
  transmute(
    patient_id = submitter_id,
    gender     = tolower(.data[["sex_at_birth"]]),
    age        = round(age_at_diagnosis / 365.25, 1),   # age_at_diagnosis is in days
    stage      = ajcc_pathologic_stage,
    T_stage    = ajcc_pathologic_t,
    N_stage    = ajcc_pathologic_n,
    M_stage    = ajcc_pathologic_m,
    race, vital_status, days_to_death, days_to_last_follow_up
  ) %>%
  distinct(patient_id, .keep_all = TRUE)
dim(clin_core)
```

```
## [1] 1098   11
```

``` r
table(clin_core$gender, useNA = "ifany")
```

```
## 
## female   male   <NA> 
##   1085     12      1
```

``` r
stopifnot(!anyDuplicated(clin_core$patient_id))
```

## 3. Clinical data B: BCR Biotab (ER / PR / HER2 by IHC)

ER/PR/HER2 IHC status is in the BCR Biotab supplement, not in the indexed clinical table.


``` r
query_clin <- GDCquery(
  project       = "TCGA-BRCA",
  data.category = "Clinical",
  data.type     = "Clinical Supplement",
  data.format   = "BCR Biotab"
)
GDCdownload(query_clin, directory = file.path(raw_dir, "GDCdata"))
biotab <- GDCprepare(query_clin, directory = file.path(raw_dir, "GDCdata"))
saveRDS(biotab, f_biotab)                         # checkpoint
```


``` r
biotab <- readRDS(f_biotab)

receptor <- biotab$clinical_patient_brca %>%
  dplyr::slice(-1, -2) %>%                                # first two rows are CDE header descriptions
  transmute(
    patient_id = bcr_patient_barcode,
    ER_status  = er_status_by_ihc,
    PR_status  = pr_status_by_ihc,
    HER2_IHC   = her2_status_by_ihc,
    histology  = histological_type
  ) %>%
  distinct(patient_id, .keep_all = TRUE)
dim(receptor)
```

```
## [1] 1097    5
```

``` r
table(receptor$ER_status, useNA = "ifany")
```

```
## 
## [Not Evaluated]   Indeterminate        Negative        Positive 
##              49               2             238             808
```

## 4. Sample selection log

Each step records number of samples / patients and the male/female split.


``` r
log_step <- function(barcodes, step, log = NULL) {
  pid <- unique(substr(barcodes, 1, 12))
  sex <- clin_core$gender[match(pid, clin_core$patient_id)]
  rbind(log, data.table(step = step, n_samples = length(barcodes), n_patients = length(pid),
                        male = sum(sex == "male", na.rm = TRUE),
                        female = sum(sex == "female", na.rm = TRUE)))
}

# Step 0: all RNA-seq aliquots in TCGA-BRCA
bc0  <- colnames(se)
slog <- log_step(bc0, "0. All TCGA-BRCA RNA-seq aliquots (GDC STAR-Counts)")

# Step 1: primary tumors only (sample type 01); drops solid tissue normal (11) and metastatic (06)
table(substr(bc0, 14, 15))
```

```
## 
##   01   06   11 
## 1111    7  113
```

``` r
bc1  <- bc0[substr(bc0, 14, 15) == "01"]
slog <- log_step(bc1, "1. Primary tumor samples (sample type 01)", slog)

# Step 2: one aliquot per patient (keep first occurrence; TCGA replicate aliquots)
bc2  <- bc1[!duplicated(substr(bc1, 1, 12))]
slog <- log_step(bc2, "2. One aliquot per patient", slog)

# Step 3: has indexed clinical data with sex
pid3 <- substr(bc2, 1, 12)
bc3  <- bc2[pid3 %in% clin_core$patient_id[!is.na(clin_core$gender)]]
slog <- log_step(bc3, "3. Clinical data with sex available", slog)

# Step 4: ER-positive by IHC (study cohort; all TCGA MBC cases are ER+)
er_pos <- receptor$patient_id[receptor$ER_status == "Positive" & !is.na(receptor$ER_status)]
bc4  <- bc3[substr(bc3, 1, 12) %in% er_pos]
slog <- log_step(bc4, "4. ER-positive by IHC (BCR Biotab)", slog)

knitr::kable(slog, caption = "Sample selection log, TCGA-BRCA")
```



Table: Sample selection log, TCGA-BRCA

|step                                                | n_samples| n_patients| male| female|
|:---------------------------------------------------|---------:|----------:|----:|------:|
|0. All TCGA-BRCA RNA-seq aliquots (GDC STAR-Counts) |      1231|       1095|   12|   1082|
|1. Primary tumor samples (sample type 01)           |      1111|       1095|   12|   1082|
|2. One aliquot per patient                          |      1095|       1095|   12|   1082|
|3. Clinical data with sex available                 |      1094|       1094|   12|   1082|
|4. ER-positive by IHC (BCR Biotab)                  |       807|        807|   12|    795|

``` r
fwrite(slog, file.path(out_dir, "TCGA_BRCA_sample_selection_log.tsv"), sep = "\t")

# Manuscript numbers: 12 male / 784 female ER+. GDC releases change slightly over time,
# so male is a hard check and female is reported for comparison.
stopifnot(slog[.N, male] == 12)
if (slog[.N, female] != 784)
  message("NOTE: female ER+ n = ", slog[.N, female], " (manuscript: 784). Document in Methods.")

# What was excluded at the ER step (for the response letter)
excl <- data.table(patient_id = substr(setdiff(bc3, bc4), 1, 12))
excl <- merge(excl, receptor[, c("patient_id", "ER_status")], by = "patient_id", all.x = TRUE)
excl <- merge(excl, clin_core[, c("patient_id", "gender")], by = "patient_id", all.x = TRUE)
excl[, .N, by = .(ER_status, gender)][order(gender, ER_status)]
```

```
##          ER_status gender     N
##             <char> <char> <int>
## 1:   Indeterminate female     2
## 2:        Negative female   237
## 3: [Not Evaluated] female    48
```

## 5. Expression matrices (all primary tumors, then ER+ cohort)


``` r
se_tumor <- se[, bc3]
stopifnot(!anyDuplicated(substr(colnames(se_tumor), 1, 12)))
colnames(se_tumor) <- substr(colnames(se_tumor), 1, 12)   # patient IDs

counts <- assay(se_tumor, "unstranded")
tpm    <- assay(se_tumor, "tpm_unstrand")

gene_info <- as.data.frame(rowData(se_tumor)) %>%
  select(gene_id, gene_name, gene_type)
gene_info$ensembl_id <- sub("\\.\\d+$", "", gene_info$gene_id)

# Symbol-level matrices: per symbol keep the Ensembl row with highest mean expression
collapse_to_symbol <- function(mat, symbols) {
  keep <- !is.na(symbols) & symbols != ""
  mat  <- mat[keep, ]; symbols <- symbols[keep]
  o    <- order(-rowMeans(mat))
  mat  <- mat[o, ]; symbols <- symbols[o]
  mat  <- mat[!duplicated(symbols), ]
  rownames(mat) <- symbols[!duplicated(symbols)]
  mat
}
counts_sym <- collapse_to_symbol(counts, gene_info$gene_name)
tpm_sym    <- collapse_to_symbol(tpm,    gene_info$gene_name)
dim(counts_sym); counts_sym[1:5, 1:3]
```

```
## [1] 59427  1094
```

```
##         TCGA-A8-A09E TCGA-A7-A0DC TCGA-A8-A08F
## MT-CO1        424249        20639       337330
## MT-ND4        329110         5320       440533
## MT-RNR2       208152      3124832       262191
## COL1A1        173022        20347       531457
## COL1A2        109465        32445       320403
```

## 6. Merge clinical tables and align with expression


``` r
clinical_full <- clin_core %>%
  left_join(receptor, by = "patient_id") %>%
  filter(patient_id %in% colnames(counts_sym)) %>%
  arrange(match(patient_id, colnames(counts_sym)))
dim(clinical_full)
```

```
## [1] 1094   15
```

``` r
stopifnot(identical(colnames(counts_sym), clinical_full$patient_id),
          identical(colnames(tpm_sym),    clinical_full$patient_id))

table(clinical_full$gender, clinical_full$ER_status, useNA = "ifany")
```

```
##         
##          [Not Evaluated] Indeterminate Negative Positive
##   female              48             2      237      795
##   male                 0             0        0       12
```

## 7. ER-positive cohort


``` r
clinical_ER <- clinical_full %>% filter(ER_status == "Positive")
counts_ER   <- counts_sym[, clinical_ER$patient_id]
tpm_ER      <- tpm_sym[,   clinical_ER$patient_id]
counts_ER_ensembl <- counts[, clinical_ER$patient_id]
tpm_ER_ensembl    <- tpm[,    clinical_ER$patient_id]

table(clinical_ER$gender)
```

```
## 
## female   male 
##    795     12
```

``` r
stopifnot(sum(clinical_ER$gender == "male") == 12)
```

## 8. Write outputs


``` r
wt <- function(x, f, rn = "GeneSymbol")
  write_tsv(as.data.frame(x) %>% rownames_to_column(rn), file.path(out_dir, f))

# ER+ cohort (study cohort)
wt(counts_ER,          "output_TCGA_BRCA_ERpos_counts.tsv")
wt(tpm_ER,             "output_TCGA_BRCA_ERpos_TPM.tsv")
wt(counts_ER_ensembl,  "output_TCGA_BRCA_ERpos_counts_ensembl.tsv", "gene_id")
wt(tpm_ER_ensembl,     "output_TCGA_BRCA_ERpos_TPM_ensembl.tsv",    "gene_id")
write_tsv(clinical_ER, file.path(out_dir, "output_TCGA_BRCA_ERpos_clinical.tsv"))

# All primary tumors (for reference / alternative cohorts)
write_tsv(clinical_full, file.path(out_dir, "output_TCGA_BRCA_all_primary_clinical.tsv"))
write_tsv(gene_info,     file.path(out_dir, "output_TCGA_BRCA_gene_annotation.tsv"))

saveRDS(list(counts_sym = counts_sym, tpm_sym = tpm_sym,
             counts_ensembl = counts, tpm_ensembl = tpm,
             clinical = clinical_full, gene_info = gene_info, selection_log = slog,
             version = paste("GDC TCGA-BRCA STAR-Counts, downloaded", Sys.Date())),
        file.path(out_dir, "TCGA_BRCA_primary.rds"))

list.files(out_dir)
```

```
## [1] "output_TCGA_BRCA_all_primary_clinical.tsv"
## [2] "output_TCGA_BRCA_ERpos_clinical.tsv"      
## [3] "output_TCGA_BRCA_ERpos_counts_ensembl.tsv"
## [4] "output_TCGA_BRCA_ERpos_counts.tsv"        
## [5] "output_TCGA_BRCA_ERpos_TPM_ensembl.tsv"   
## [6] "output_TCGA_BRCA_ERpos_TPM.tsv"           
## [7] "output_TCGA_BRCA_gene_annotation.tsv"     
## [8] "TCGA_BRCA_primary.rds"                    
## [9] "TCGA_BRCA_sample_selection_log.tsv"
```


``` r
sessionInfo()
```

```
## R version 4.5.3 (2026-03-11)
## Platform: aarch64-apple-darwin20
## Running under: macOS Tahoe 26.6.2
## 
## Matrix products: default
## BLAS:   /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRblas.0.dylib 
## LAPACK: /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## time zone: America/Los_Angeles
## tzcode source: internal
## 
## attached base packages:
## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
## [8] base     
## 
## other attached packages:
##  [1] lubridate_1.9.5             forcats_1.0.1              
##  [3] stringr_1.6.0               dplyr_1.2.1                
##  [5] purrr_1.2.2                 readr_2.2.0                
##  [7] tidyr_1.3.2                 tibble_3.3.1               
##  [9] ggplot2_4.0.3               tidyverse_2.0.0            
## [11] data.table_1.18.6.1         SummarizedExperiment_1.40.0
## [13] Biobase_2.70.0              GenomicRanges_1.62.1       
## [15] Seqinfo_1.0.0               IRanges_2.44.0             
## [17] S4Vectors_0.48.1            BiocGenerics_0.56.0        
## [19] generics_0.1.4              MatrixGenerics_1.22.0      
## [21] matrixStats_1.5.0           TCGAbiolinks_2.38.0        
## 
## loaded via a namespace (and not attached):
##  [1] tidyselect_1.2.1            farver_2.1.2               
##  [3] blob_1.3.0                  filelock_1.0.3             
##  [5] Biostrings_2.78.0           S7_0.2.2                   
##  [7] fastmap_1.2.0               BiocFileCache_3.0.0        
##  [9] XML_3.99-0.24               digest_0.6.39              
## [11] timechange_0.4.0            lifecycle_1.0.5            
## [13] KEGGREST_1.50.0             RSQLite_3.53.3             
## [15] magrittr_2.0.5              compiler_4.5.3             
## [17] rlang_1.3.0                 sass_0.4.10                
## [19] progress_1.2.3              tools_4.5.3                
## [21] yaml_2.3.12                 knitr_1.51                 
## [23] prettyunits_1.2.0           S4Arrays_1.10.1            
## [25] bit_4.6.0                   curl_8.0.0                 
## [27] DelayedArray_0.36.1         plyr_1.8.9                 
## [29] xml2_1.6.0                  RColorBrewer_1.1-3         
## [31] abind_1.4-8                 withr_3.0.3                
## [33] grid_4.5.3                  scales_1.4.0               
## [35] biomaRt_2.66.2              cli_3.6.6                  
## [37] rmarkdown_2.32              crayon_1.5.3               
## [39] otel_0.2.0                  rstudioapi_0.19.0          
## [41] httr_1.4.9                  tzdb_0.5.0                 
## [43] DBI_1.3.0                   cachem_1.1.0               
## [45] parallel_4.5.3              rvest_1.0.5                
## [47] AnnotationDbi_1.72.0        TCGAbiolinksGUI.data_1.30.0
## [49] XVector_0.50.0              vctrs_0.7.3                
## [51] Matrix_1.7-6                jsonlite_2.0.0             
## [53] hms_1.1.4                   bit64_4.8.6                
## [55] jquerylib_0.1.4             glue_1.8.1                 
## [57] stringi_1.8.9               gtable_0.3.6               
## [59] pillar_1.11.1               htmltools_0.5.9            
## [61] R6_2.6.1                    dbplyr_2.6.0               
## [63] httr2_1.3.0                 vroom_1.7.1                
## [65] evaluate_1.0.5              lattice_0.23-1             
## [67] png_0.1-9                   memoise_2.0.1              
## [69] bslib_0.12.0                Rcpp_1.1.2                 
## [71] SparseArray_1.10.10         downloader_0.4.1           
## [73] xfun_0.60                   pkgconfig_2.0.3
```
