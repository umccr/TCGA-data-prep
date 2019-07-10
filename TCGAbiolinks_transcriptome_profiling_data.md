# Processing harmonized transcriptome profiling data

Pipeline for downloading HARMONIZED transcriptome profiling data from [TCGA](https://cancergenome.nih.gov) and [TARGET](https://ocg.cancer.gov/programs/target/research) projects using *[TCGAbiolinks](https://bioconductor.org/packages/release/bioc/vignettes/TCGAbiolinks/inst/doc/index.html)* R package. It uses the [GDC Data Transfer Tool Client](https://gdc.cancer.gov/access-data/gdc-data-transfer-tool) to download the data from the [Genomic Data Commons](https://gdc.cancer.gov/about-gdc) (GDC) portal. The script outputs normalised expression (FPKM or FPKM-UQ) or raw count matrix for RNA-seq data for user-defined tissue types along with associated clinical information.

The TCGA genomic data harmonization is [here](https://gdc.cancer.gov/about-data/data-harmonization-and-generation/genomic-data-harmonization-0#Overview) and the mRNA analysis pipeline is described [here](https://gdc-docs.nci.nih.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline).
<br />
<br />

## Table of contents

<!-- vim-markdown-toc GFM -->
* [Installation](#installation)
* [Arguments](#arguments)
* [Example of use](#example-of-use)
  * [Command line](#command-line)
  * [Output data directory structure](#output-data-directory-structure)
  * [Files description](#files-description)
* [Arguments options](#arguments-options)
  * [--out_dir](#--out_dir)
  * [--project_id](#--project_id)
  * [--tissue](#--tissue)
  * [--workflow](#--workflow)


<!-- vim-markdown-toc -->
<br>

## Installation

* **Stable version** of [TCGAbiolinks](http://bioconductor.org/packages/release/bioc/vignettes/TCGAbiolinks/inst/doc/index.html) R package from [Bioconductor](http://bioconductor.org/packages/release/bioc/html/TCGAbiolinks.html)

```
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("TCGAbiolinks")
```

* NOTE, in case of any issues with the stable version, one can try using the **development version** from [GitHub](https://github.com/BioinformaticsFMRP/TCGAbiolinks)

```
devtools::install_github('BioinformaticsFMRP/TCGAbiolinks')
```

>For instance, that may be the case when [biomaRt](https://bioconductor.org/packages/release/bioc/html/biomaRt.html) is updated and significant changes are introduced to *Ensembl BioMart* database (see [this GitHub post](https://github.com/BioinformaticsFMRP/TCGAbiolinks/issues/250)).


## Arguments

Argument | Short | Description
------------ | ------------ | ------------
--out_dir | -o| Directory to which the data is to be downloaded
--project_id | -p| ID of the TCGA/TARGET project to download
--tissue | -t| Tissue types to be considered for download
--workflow | -w| Workflow from which the data is to be downloaded
<br />

## Example of use

Use the [TCGAbiolinks_transcriptome_profiling_data.R](https://github.com/umccr/TCGA-data-prep/blob/master/TCGAbiolinks_transcriptome_profiling_data.R) script to download read count matrix for pancreatic cancer (PAAD) project.

### Command line

```
Rscript TCGAbiolinks_transcriptome_profiling_data.R --out_dir TCGA/PAAD --project_id TCGA-PAAD --tissue 1,11 --workflow Counts
```

<br />

### Output data directory structure

```
TCGA
|
|____PAAD
  |
  |____transcriptome_profiling
    |
    |____Counts
      |____Counts.exp
      |____Counts_boxplot.pdf
      |____Counts_clinical_info.txt
      |____Counts_samples.txt
      |____gdc-client
      |____gdc-client_v1.1.0_OSX_x64.zip
      |____gdc_manifest.txt
      |____R_parameters.txt
      |____GDCdata
        |
        |____TCGA-PAAD
          |
          |____harmonized
            |
            |____Transcriptome_Profiling
            | |
            | |____Gene_Expression_Quantification
            |   |____…
            |   |____…
            |
            |____Clinical
              |
              |____Clinical_Supplement
                |____…
                |____…
```
<br />

### Files description

File | Description
------------ | ------------
Counts.exp | Read count data matrix
Counts_boxplot.pdf | Box plot of read counts per sample
Counts_clinical_info.txt | Samples and associated clinical annotation
R_parameters.txt | User-defined parameters used for the script execution
Gene_Expression_Quantification | Folder with compressed 'txt' files containing expression values for each sample
Clinical_Supplement | Folder with 'xml' files including clinical information for each sample
<br />


## Arguments options

### --out_dir

Local workspace. This is the directory to which the data will be downloaded and stored.

<br />

### --project_id

Available TCGA/TARGET project IDs are:

Project ID | Name
------------ | ------------
TCGA-SARC | Sarcoma
TCGA-MESO | Mesothelioma
TCGA-READ | Rectum Adenocarcinoma
TCGA-KIRP | Kidney Renal Papillary Cell Carcinoma
TARGET-NBL | Neuroblastoma
TCGA-PAAD  | Pancreatic Adenocarcinoma
TCGA-GBM   | Glioblastoma Multiforme
TCGA-ACC   | Adrenocortical Carcinoma
TARGET-OS  | Osteosarcoma
TCGA-CESC  | Cervical Squamous Cell Carcinoma and Endocervical Adenocarcinoma
TARGET-RT  | Rhabdoid Tumour
TCGA-BRCA  | Breast Invasive Carcinoma
TCGA-ESCA  | Esophageal Carcinoma
TCGA-DLBC  | Lymphoid Neoplasm Diffuse Large B-cell Lymphoma
TCGA-KICH  | Kidney Chromophobe
TCGA-KIRC  | Kidney Renal Clear Cell Carcinoma
TCGA-UVM   | Uveal Melanoma
TARGET-AML | Acute Myeloid Leukaemia
TCGA-LAML  | Acute Myeloid Leukaemia
TCGA-SKCM  | Skin Cutaneous Melanoma
TCGA-PCPG  | Pheochromocytoma and Paraganglioma
TCGA-COAD  | Colon Adenocarcinoma
TCGA-UCS   | Uterine Carcinosarcoma
TCGA-LUSC  | Lung Squamous Cell Carcinoma
TCGA-LGG   | Brain Lower Grade Glioma
TCGA-HNSC  | Head and Neck Squamous Cell Carcinoma
TCGA-TGCT  | Testicular Germ Cell Tumours
TARGET-CCSK  | Clear Cell Sarcoma of the Kidney
TCGA-THCA  | Thyroid Carcinoma
TCGA-LIHC  | Liver Hepatocellular Carcinoma
TCGA-BLCA  | Bladder Urothelial Carcinoma
TCGA-UCEC  | Uterine Corpus Endometrial Carcinoma
TARGET-WT  | High-Risk Wilms Tumour
TCGA-PRAD  | Prostate Adenocarcinoma
TCGA-OV    | Ovarian Serous Cystadenocarcinoma
TCGA-THYM  | Thymoma
TCGA-CHOL  | Cholangiocarcinoma
TCGA-STAD  | Stomach Adenocarcinoma
TCGA-LUAD  | Lung Adenocarcinoma
<br />

### --tissue

Multiple tissue types are allowed. Each tissue type is expected to be separated by comma. Type '*all*' for all listed tissue types to be considered for download. Available options are:

Tissue code | Letter code | Definition
------------ | ------------ | ------------
1 | TP  | Primary solid Tumour
2 | TR  | Recurrent Solid Tumour
3 | TB  | Primary Blood Derived Cancer - Peripheral Blood
4 | TRBM | Recurrent Blood Derived Cancer - Bone Marrow
5 | TAP | Additional - New Primary
6 | TM | Metastatic
7 | TAM | Additional Metastatic
8 | THOC | Human Tumour Original Cells
9 | TBM | Primary Blood Derived Cancer - Bone Marrow
10 | NB | Blood Derived Normal
11 | NT | Solid Tissue Normal
12 | NBC | Buccal Cell Normal
13 | NEBV | EBV Immortalised Normal
14 | NBM | Bone Marrow Normal
20 | CELLC | Control Analyte
40 | TRB | Recurrent Blood Derived Cancer - Peripheral Blood
50 | CELL | Cell Lines
60 | XP | Primary Xenograft Tissue
61 | XCL | Cell Line Derived Xenograft Tissue
All | --- | All available tissue types
<br />

### --workflow

Data from three [workflows](https://gdc-docs.nci.nih.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline) are available:

Workflow | Definition
------------ | ------------
Counts | Raw Read Counts - the number of reads aligned to each protein-coding gene, calculated by *HT-Seq* (**default**)
FPKM | Normalised expression value that takes into account each protein-coding gene length and the number of reads mappable to all protein-coding genes
FPKM-UQ | Normalised raw read count in which gene expression values, in FPKM, are divided by the 75th percentile value
<br />

## Note
Make sure that R version >= 3.3 is installed. For older versions the *TCGAbiolinks* uses different functions starting with "TCGA" rather than "GDC" since the data were moved from DCC server to NCI Genomic Data Commons (GDC).
Make sure that the [newest *TCGAbiolinks* package](https://github.com/BioinformaticsFMRP/TCGAbiolinks) package is installed.

```
devtools::install_github(repo = "BioinformaticsFMRP/TCGAbiolinks")
```
