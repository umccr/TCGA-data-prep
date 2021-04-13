# Processing legacy gene expression data

Pipeline for downloading LEGACY (genome reference `hg19`) gene expression data from [TCGA](https://cancergenome.nih.gov) and [TARGET](https://ocg.cancer.gov/programs/target/research) projects using *[TCGAbiolinks](http://bioconductor.org/packages/release/bioc/vignettes/TCGAbiolinks/inst/doc/tcgaBiolinks.html)* R package. The script outputs normalised expression or raw count matrix for RNA-seq data for user-defined tissue types along with associated clinical information.

<br />
<br />

### Arguments
1. **Local workspace**. This is the directory to which the data will be downloaded and stored.

2. **Project ID**. Available TCGA/TARGET project IDs are:

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

3. **Tissue code**. Tissue types to be considered for download. Each tissue type is expected to be separated by comma. Type '*all*' for all listed tissue types to be considered for download. Available options are:

      Tissue code | Letter code | Definition
      ------------ | ------------ | ------------
      01 | TP  | Primary solid Tumour
      02 | TR  | Recurrent Solid Tumour
      03 | TB  | Primary Blood Derived Cancer - Peripheral Blood
      04 | TRBM | Recurrent Blood Derived Cancer - Bone Marrow
      05 | TAP | Additional - New Primary
      06 | TM | Metastatic
      07 | TAM | Additional Metastatic
      08 | THOC | Human Tumour Original Cells
      09 | TBM | Primary Blood Derived Cancer - Bone Marrow
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

4. **Normalised data**. If "*norm*" is specified, the upper quartile normalised RSEM count estimates ("[*normalized_count*](https://wiki.nci.nih.gov/display/tcga/rnaseq+version+2)") are provided. Otherwise, read count data reflecting the number of reads mapping to each gene are provided (**default**).
<br />

###   Command line use example
```
      R --file=./TCGAbiolinks_gene_expression_data.R --args "TCGA_data/TCGA-PAAD" "TCGA-PAAD" "1,11" "norm"
```
<br />

###   Output data directory structure
```
TCGA_data
|
|____TCGA-PAAD
  |
  |____gene_expression
    |
    |____normalized
      |____normalized.exp
      |____normalized_boxplot.pdf
      |____normalized_clinical_info.txt
      |____normalized_samples.txt
      |____gdc-client
      |____gdc-client_v1.2.0_OSX_x64.zip
      |____gdc_manifest.txt
      |____R_parameters.txt
      |____GDCdata
        |
        |____TCGA-PAAD
          |
          |____legacy
          | |
          | |____Gene_expression
          |   |
          |   |____Gene_expression_quantification
          |     |____…
          |     |____…
          |
          |____harmonized
            |
            |____Clinical
              |
              |____Clinical_Supplement
                |____…
                |____…
```
<br />

###   Files description
File | Description
------------ | ------------
normalized.exp | Gene expression matrix
normalized_boxplot.pdf | Box plot of normalised log2 expression data per sample
normalized_clinical_info.txt | Samples and associated clinical annotation
R_parameters.txt | User-defined parameters used for the script execution
Gene_expression_quantification | Folder with compressed 'txt' files containing expression values for each sample
Clinical_Supplement | Folder with 'xml' files including clinical information for each sample
<br />

### Note
Make sure that R version >= 3.3 is installed. For older versions the *TCGAbiolinks* uses different functions starting with "TCGA" rather than "GDC" since the data were moved from DCC server to NCI Genomic Data Commons (GDC).
Make sure that the [newest *TCGAbiolinks* package](https://github.com/BioinformaticsFMRP/TCGAbiolinks) package is installed.

```
devtools::install_github(repo = "BioinformaticsFMRP/TCGAbiolinks")
```
