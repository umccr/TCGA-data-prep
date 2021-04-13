# Processing HARMONIZED DNA methylation data

## NOTE: Needs updates!!!

Pipeline for downloading HARMONIZED DNA methylation data from [TCGA](https://cancergenome.nih.gov) and [TARGET](https://ocg.cancer.gov/programs/target/research) projects using *[TCGAbiolinks](http://bioconductor.org/packages/release/bioc/vignettes/TCGAbiolinks/inst/doc/tcgaBiolinks.html)* R package. The script outputs matrix with DNA methylation data for user-defined tissue types along with associated clinical information.

The TCGA genomic data harmonization is [here](https://gdc.cancer.gov/about-data/data-harmonization-and-generation/genomic-data-harmonization-0#Overview) and the the DNA methylation liftover pipeline is described [here](https://gdc-docs.nci.nih.gov/Data/Bioinformatics_Pipelines/Methylation_LO_Pipeline).

###### Note

Genome of reference `hg38` is used as default.

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

4. **Platform**. Available data may have been produced using one of two platforms:

      Platform | Full name
      ------------ | ------------
      HM27 | Illumina Infinium Human Methylation 27
      HM450 | Illumina Infinium Human HumanMethylation 450
<br />

###   Command line use example
```
      R --file=./TCGAbiolinks_DNA_methylation_data.R --args "TCGA_data/TCGA-PAAD" "TCGA-PAAD" "1,11" "HM450"
```
<br />

###   Output data directory structure
```
TCGA_data
|
|____TCGA-PAAD
  |
  |____DNA_methylation
    |
    |____HM450
      |____HM450.txt
      |____HM450_boxplot.pdf
      |____HM450_clinical_info.txt
      |____HM450_samples.txt
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
            |____DNA_Methylation
            | |
            | |____Methylation_Beta_Value
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

###   Files description
File | Description
------------ | ------------
HM450.txt | DNA methylation (Beta values) matrix
HM450_boxplot.pdf | Box plot of DNA methylation data (Beta values) per sample
HM450_clinical_info.txt | Samples and associated clinical annotation
R_parameters.txt | User-defined parameters used for the script execution
DNA_Methylation | Folder with compressed 'txt' files containing DNA methylation (Beta) values for each sample
Clinical_Supplement | Folder with 'xml' files including clinical information for each sample
<br />

### Note
Make sure that R version >= 3.3 is installed. For older versions the *TCGAbiolinks* uses different functions starting with "TCGA" rather than "GDC" since the data were moved from DCC server to NCI Genomic Data Commons (GDC).
Make sure that the [newest *TCGAbiolinks* package](https://github.com/BioinformaticsFMRP/TCGAbiolinks) is installed.

```
devtools::install_github(repo = "BioinformaticsFMRP/TCGAbiolinks")
```
