# TCGA data preparation for analyses

Pipelines for downloading and preparing different types of [TCGA](https://cancergenome.nih.gov) and [TARGET](https://ocg.cancer.gov/programs/target/research) data for analyses using the *[TCGAbiolinks](https://bioconductor.org/packages/release/bioc/vignettes/TCGAbiolinks/inst/doc/index.html)* R package:

- [TCGAbiolinks_transcriptome_profiling_data](TCGAbiolinks_transcriptome_profiling_data.md)
pipeline for processing HARMONIZED transcriptome profiling data

- [TCGAbiolinks_gene_expression_data](TCGAbiolinks_gene_expression_data.md)
pipeline for processing LEGACY gene expression data

- [TCGAbiolinks_DNA_methylation_data](TCGAbiolinks_DNA_methylation_data.md)
pipeline for processing HARMONIZED DNA methylation data

<br />

The [TCGAbiolinks](http://bioconductor.org/packages/release/bioc/vignettes/TCGAbiolinks/inst/doc/index.html) package accesses the [National Cancer Institute](https://www.cancer.gov/) (NCI) [GDC](https://gdc.cancer.gov/about-data/publications/pancanatlas) thorough its [GDC Application Programming Interface](https://gdc.cancer.gov/developers/gdc-application-programming-interface-api) (API) to download the most recent data files.