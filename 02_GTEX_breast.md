---
title: "Generate gene-reads / TPM matrices for GTEx v8 male and female breast"
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
library(data.table)
knitr::opts_knit$set(root.dir = ".")

raw_dir <- "Raw_Data"
out_dir <- "Output"
dir.create(out_dir, showWarnings = FALSE)

f_attr  <- file.path(raw_dir, "GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")
f_pheno <- file.path(raw_dir, "GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt")
f_reads <- file.path(raw_dir, "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz")
f_tpm   <- file.path(raw_dir, "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz")
stopifnot(file.exists(f_attr, f_pheno, f_reads, f_tpm))
```

## 1. Sample attributes: all breast aliquots


``` r
data1 <- fread(f_attr, quote = "")
dim(data1)
```

```
## [1] 22951    63
```

``` r
breast_sample <- data1[SMTS == "Breast", SAMPID]
length(breast_sample)
```

```
## [1] 480
```

``` r
stopifnot(length(breast_sample) == 480)
### There are 480 breast samples (aliquots), including male and female ###
```

## 2. Subject phenotypes: donor sex


``` r
data2 <- fread(f_pheno, quote = "")
dim(data2)
```

```
## [1] 980   4
```

``` r
### There are 980 donors in total ###
male_donor   <- data2[SEX == 1, SUBJID]
female_donor <- data2[SEX == 2, SUBJID]
length(male_donor); length(female_donor)
```

```
## [1] 653
```

```
## [1] 327
```

``` r
stopifnot(length(male_donor) == 653, length(female_donor) == 327)
### 653 male donors, 327 female donors ###

# Sample ID "GTEX-XXXXX-0226-SM-XXXXX" -> donor ID "GTEX-XXXXX"
breast_donor <- sub("^(GTEX-[^-]+)-.*$", "\\1", breast_sample)
breast_sample_male   <- breast_sample[breast_donor %in% male_donor]
breast_sample_female <- breast_sample[breast_donor %in% female_donor]
length(breast_sample_male); length(breast_sample_female)
```

```
## [1] 302
```

```
## [1] 178
```

``` r
stopifnot(length(breast_sample_male) == 302, length(breast_sample_female) == 178)
### Among 480 breast aliquots: 302 male, 178 female ###
```

## 3. Read expression matrices (breast columns only)


``` r
read_gct_breast <- function(file, samples) {
  cols <- names(fread(file, skip = 2, nrows = 0))
  idx  <- which(cols %in% c("Name", "Description", samples))
  dt   <- fread(file, skip = 2, select = idx)
  m    <- as.data.frame(dt[, -"Name"])       
  rownames(m) <- dt$Name
  m
}

rna_breast <- read_gct_breast(f_reads, c(breast_sample_male, breast_sample_female))
tpm_breast <- read_gct_breast(f_tpm,   c(breast_sample_male, breast_sample_female))
dim(rna_breast); rna_breast[1:5, 1:5]
```

```
## [1] 56200   460
```

```
##                   Description GTEX-1117F-2826-SM-5GZXL GTEX-111YS-1926-SM-5GICC
## ENSG00000223972.5     DDX11L1                        0                        1
## ENSG00000227232.5      WASH7P                      286                      135
## ENSG00000278267.1   MIR6859-1                        0                        0
## ENSG00000243485.5 MIR1302-2HG                        0                        0
## ENSG00000237613.2     FAM138A                        0                        0
##                   GTEX-1122O-1226-SM-5H113 GTEX-117XS-1926-SM-5GICO
## ENSG00000223972.5                        0                        1
## ENSG00000227232.5                      110                      246
## ENSG00000278267.1                        0                        0
## ENSG00000243485.5                        0                        0
## ENSG00000237613.2                        0                        0
```

``` r
stopifnot(identical(rownames(rna_breast), rownames(tpm_breast)),
          identical(colnames(rna_breast), colnames(tpm_breast)))
### 56200 genes; only samples present in the V8 analysis freeze have a column ###
```

## 4. Split by sex and check sample numbers


``` r
in_matrix     <- colnames(rna_breast)[-1]
male_sample_with_RNA     <- breast_sample_male[breast_sample_male %in% in_matrix]
female_sample_with_RNA   <- breast_sample_female[breast_sample_female %in% in_matrix]
length(male_sample_with_RNA); length(female_sample_with_RNA)
```

```
## [1] 291
```

```
## [1] 168
```

``` r
stopifnot(length(male_sample_with_RNA) == 291, length(female_sample_with_RNA) == 168)
### 291 male and 168 female breast samples with RNA-seq data ###

# 21 aliquots dropped: all flagged SMAFRZE == "EXCLUDE" by GTEx QC (11 male, 10 female)
dropped <- data1[SAMPID %in% setdiff(c(breast_sample_male, breast_sample_female), in_matrix)]
dropped[, .N, by = SMAFRZE]
```

```
##    SMAFRZE     N
##     <char> <int>
## 1: EXCLUDE    21
```

``` r
dropped_sex <- merge(dropped[, .(SAMPID, SUBJID = sub("^(GTEX-[^-]+)-.*$", "\\1", SAMPID), SMAFRZE)],
                     data2[, .(SUBJID, SEX)], by = "SUBJID")
dropped_sex[, .N, by = .(SMAFRZE, SEX)]
```

```
##    SMAFRZE   SEX     N
##     <char> <int> <int>
## 1: EXCLUDE     2    10
## 2: EXCLUDE     1    11
```

``` r
stopifnot(setequal(in_matrix, data1[SMTS == "Breast" & SMAFRZE == "RNASEQ", SAMPID]))

rna_male_breast   <- rna_breast[, c("Description", male_sample_with_RNA)]
rna_female_breast <- rna_breast[, c("Description", female_sample_with_RNA)]
tpm_male_breast   <- tpm_breast[, c("Description", male_sample_with_RNA)]
tpm_female_breast <- tpm_breast[, c("Description", female_sample_with_RNA)]
```

## 5. Sample metadata


``` r
meta <- merge(
  data1[SAMPID %in% in_matrix, .(SAMPID, SUBJID = sub("^(GTEX-[^-]+)-.*$", "\\1", SAMPID),
                                  SMTSD, SMRIN, SMTSISCH, SMNABTCH, SMGEBTCH, SMAFRZE)],
  data2[, .(SUBJID, SEX, AGE, DTHHRDY)], by = "SUBJID")
meta[, SEX_label := ifelse(SEX == 1, "Male", "Female")]
meta <- meta[match(in_matrix, SAMPID)]
stopifnot(!anyDuplicated(meta$SUBJID))   # one breast sample per donor
meta[, .N, by = SEX_label]
```

```
##    SEX_label     N
##       <char> <int>
## 1:    Female   168
## 2:      Male   291
```

## 6. Write outputs


``` r
wt <- function(x, f) write.table(x, file.path(out_dir, f), sep = "\t",
                                 row.names = TRUE, col.names = NA, quote = FALSE)

# Same files/format as 2020
wt(rna_male_breast,   "output_GTEx_RNA_male_breast.txt")
wt(rna_female_breast, "output_GTEx_RNA_female_breast.txt")

# New: TPM, combined, metadata
wt(tpm_male_breast,   "output_GTEx_TPM_male_breast.txt")
wt(tpm_female_breast, "output_GTEx_TPM_female_breast.txt")
wt(rna_breast,        "output_GTEx_RNA_breast_all.txt")
wt(tpm_breast,        "output_GTEx_TPM_breast_all.txt")
fwrite(meta, file.path(out_dir, "output_GTEx_breast_metadata.txt"), sep = "\t")

saveRDS(list(counts = rna_breast, tpm = tpm_breast, meta = meta,
             version = "GTEx v8, RNASeQCv1.1.9, GENCODE v26"),
        file.path(out_dir, "GTEx_v8_breast.rds"))

list.files(out_dir)
```

```
## [1] "GTEx_v8_breast.rds"                "output_GTEx_breast_metadata.txt"  
## [3] "output_GTEx_RNA_breast_all.txt"    "output_GTEx_RNA_female_breast.txt"
## [5] "output_GTEx_RNA_male_breast.txt"   "output_GTEx_TPM_breast_all.txt"   
## [7] "output_GTEx_TPM_female_breast.txt" "output_GTEx_TPM_male_breast.txt"
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
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] data.table_1.18.6.1
## 
## loaded via a namespace (and not attached):
##  [1] digest_0.6.39     R6_2.6.1          fastmap_1.2.0     xfun_0.60        
##  [5] cachem_1.1.0      R.utils_2.13.0    knitr_1.51        htmltools_0.5.9  
##  [9] rmarkdown_2.32    lifecycle_1.0.5   cli_3.6.6         R.methodsS3_1.8.2
## [13] sass_0.4.10       jquerylib_0.1.4   compiler_4.5.3    R.oo_1.27.1      
## [17] rstudioapi_0.19.0 tools_4.5.3       evaluate_1.0.5    bslib_0.12.0     
## [21] yaml_2.3.12       otel_0.2.0        jsonlite_2.0.0    rlang_1.3.0
```
