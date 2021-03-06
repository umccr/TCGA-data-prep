################################################################################
#
#   File name: TCGAbiolinks_DNA_methylation_data.R
#
#   Authors: Jacek Marzec ( j.marzec@qmul.ac.uk )
#
#   Barts Cancer Institute
#   Queen Mary, University of London
#   Charterhouse Square, London EC1M 6BQ
#
################################################################################

################################################################################
#
#   Description: Pipeline for downloading HARMONIZED DNA methylation DATA from TCGA ( https://cancergenome.nih.gov/ ) and TARGET ( https://ocg.cancer.gov/programs/target/research ) projects using TCGAbiolinks package ( http://bioconductor.org/packages/release/bioc/vignettes/TCGAbiolinks/inst/doc/tcgaBiolinks.html ). The script outputs matrix with DNA methylation data for user-defined tissue types along with associated clinical information. The TCGA genomic data harmonization is here: https://gdc.cancer.gov/about-data/data-harmonization-and-generation/genomic-data-harmonization-0#Overview ; and the DNA methylation liftover pipeline is described here: https://gdc-docs.nci.nih.gov/Data/Bioinformatics_Pipelines/Methylation_LO_Pipeline/
#
#   Note: Make sure that R version >= 3.3 is installed. For older versions the TCGAbiolinks uses different functions starting with "TCGA" rather than "GDC" since the data were moved from DCC server to NCI Genomic Data Commons (GDC).
#   Make sure that the newest TCGAbiolinks package is installed ( https://github.com/BioinformaticsFMRP/TCGAbiolinks )
#
#   devtools::install_github(repo = "BioinformaticsFMRP/TCGAbiolinks")
#
#   Arguments:
#
#============> [1] Local workspace. This is the directory to which the data will be downloaded and stored
#
#
#============> [2] Project ID. Available TCGA/TARGET project IDs are:
#
#               Project ID      (Name)
#
#               TCGA-SARC       (Sarcoma)
#               TCGA-MESO       (Mesothelioma)
#               TCGA-READ       (Rectum Adenocarcinoma)
#               TCGA-KIRP       (Kidney Renal Papillary Cell Carcinoma)
#               TARGET-NBL      (Neuroblastoma)
#               TCGA-PAAD       (Pancreatic Adenocarcinoma)
#               TCGA-GBM        (Glioblastoma Multiforme)
#               TCGA-ACC        (Adrenocortical Carcinoma)
#               TARGET-OS       (Osteosarcoma)
#               TCGA-CESC       (Cervical Squamous Cell Carcinoma and Endocervical Adenocarcinoma)
#               TARGET-RT       (Rhabdoid Tumor)
#               TCGA-BRCA       (Breast Invasive Carcinoma)
#               TCGA-ESCA       (Esophageal Carcinoma)
#               TCGA-DLBC       (Lymphoid Neoplasm Diffuse Large B-cell Lymphoma)
#               TCGA-KICH       (Kidney Chromophobe)
#               TCGA-KIRC       (Kidney Renal Clear Cell Carcinoma)
#               TCGA-UVM        (Uveal Melanoma)
#               TARGET-AML      (Acute Myeloid Leukemia)
#               TCGA-LAML       (Acute Myeloid Leukemia)
#               TCGA-SKCM       (Skin Cutaneous Melanoma)
#               TCGA-PCPG       (Pheochromocytoma and Paraganglioma)
#               TCGA-COAD       (Colon Adenocarcinoma)
#               TCGA-UCS        (Uterine Carcinosarcoma)
#               TCGA-LUSC       (Lung Squamous Cell Carcinoma)
#               TCGA-LGG        (Brain Lower Grade Glioma)
#               TCGA-HNSC       (Head and Neck Squamous Cell Carcinoma)
#               TCGA-TGCT       (Testicular Germ Cell Tumors)
#               TARGET-CCSK     (Clear Cell Sarcoma of the Kidney)
#               TCGA-THCA       (Thyroid Carcinoma)
#               TCGA-LIHC       (Liver Hepatocellular Carcinoma)
#               TCGA-BLCA       (Bladder Urothelial Carcinoma)
#               TCGA-UCEC       (Uterine Corpus Endometrial Carcinoma)
#               TARGET-WT       (High-Risk Wilms Tumor)
#               TCGA-PRAD       (Prostate Adenocarcinoma)
#               TCGA-OV         (Ovarian Serous Cystadenocarcinoma)
#               TCGA-THYM       (Thymoma)
#               TCGA-CHOL       (Cholangiocarcinoma)
#               TCGA-STAD       (Stomach Adenocarcinoma)
#               TCGA-LUAD       (Lung Adenocarcinoma)
#
#
#============> [3] Tissue code. Tissue types to be considered for download. Each tissue type is expected to be separated by comma. Type 'all' for all listed tissue types to be considered for download. Available options are:
#
#               Tissue   (Letter    (Definition)
#               code      code)
#
#               1        (TP)       (Primary solid Tumor)
#               2        (TR)       (Recurrent Solid Tumor)
#               3        (TB)       (Primary Blood Derived Cancer - Peripheral Blood)
#               4        (TRBM)     (Recurrent Blood Derived Cancer - Bone Marrow)
#               5        (TAP)      (Additional - New Primary)
#               6        (TM)       (Metastatic)
#               7        (TAM)      (Additional Metastatic)
#               8        (THOC)     (Human Tumor Original Cells)
#               9        (TBM)      (Primary Blood Derived Cancer - Bone Marrow)
#               10       (NB)       (Blood Derived Normal)
#               11       (NT)       (Solid Tissue Normal)
#               12       (NBC)      (Buccal Cell Normal)
#               13       (NEBV)     (EBV Immortalized Normal)
#               14       (NBM)      (Bone Marrow Normal)
#               20       (CELLC)    (Control Analyte)
#               40       (TRB)      (Recurrent Blood Derived Cancer - Peripheral Blood)
#               50       (CELL)     (Cell Lines)
#               60       (XP)       (Primary Xenograft Tissue)
#               61       (XCL)      (Cell Line Derived Xenograft Tissue)
#               All       ---       (All available tissue types)
#
#
#============> [4] Platform. Available data may have been produced using one of two platforms:
#
#               Platform  (Full name)
#
#               HM27      (Illumina Infinium Human Methylation 27)
#               HM450     (Illumina Infinium Human HumanMethylation 450; used as DEFAULT)
#
#
#
#   Command line use example: R --file=./TCGAbiolinks_DNA_methylation_data.R --args "TCGA_data/TCGA-PAAD" "TCGA-PAAD" "1,11" "HM450"
#
#
################################################################################

##### Clear workspace
rm(list=ls())
##### Close any open graphics devices
graphics.off()

#===============================================================================
#    Functions
#===============================================================================

##### Create 'not in' operator
"%!in%" <- function(x,table) match(x,table, nomatch = 0) == 0


##### Prepare object to write into a file
prepare2write <- function (x) {

	x2write <- cbind(rownames(x), x)
    colnames(x2write) <- c("",colnames(x))
	return(x2write)
}


#===============================================================================
#    Load libraries
#===============================================================================

library("TCGAbiolinks")
library("SummarizedExperiment")


#===============================================================================
#    Main
#===============================================================================

args <- commandArgs()

ProjectDir <- args[4]
ProjectID <- args[5]
Tissue <- args[6]
Tissue <- gsub("\\s","", Tissue)
Tissue <- sort(unlist(strsplit(Tissue, split=',', fixed=TRUE)))

if ( Tissue == "All" || Tissue == "all" ) {

    Tissue <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,20,40,50,60,61)
}


##### Define the Infinium Human Methylation platfrom
if ( is.na(args[7]) ) {

    cat( paste("Platform type information is missing\n", sep=" ") )
    q()

} else if ( args[7] == "HM27" || args[7] == "hm27" ) {

    platform <- "Illumina Human Methylation 27"

    ##### Make the results directory more specific
    ProjectDir <- paste(ProjectDir, "DNA_methylation/HM27", sep="/" )

} else if ( args[7] == "HM450" || args[7] == "hm450" ) {

    platform <- "Illumina Human Methylation 450"

    ##### Make the results directory more specific
    ProjectDir <- paste(ProjectDir, "DNA_methylation/HM450", sep="/" )

} else {

    cat( paste("The platform type", args[7], "is not supported\n", sep=" ") )
    q()
}

##### Get the local project name
ProjectName <-  unlist(strsplit(ProjectDir, split='/', fixed=TRUE))
ProjectName <-  ProjectName[length(ProjectName)]

##### Store available project IDs in a vector
ProjectIDs <- c("TCGA-SARC","TCGA-MESO","TCGA-READ","TCGA-KIRP","TARGET-NBL","TCGA-PAAD","TCGA-GBM","TCGA-ACC","TARGET-OS","TCGA-CESC","TARGET-RT","TCGA-BRCA","TCGA-ESCA","TCGA-DLBC","TCGA-KICH","TCGA-KIRC","TCGA-UVM","TARGET-AML","TCGA-LAML","TCGA-SKCM","TCGA-PCPG","TCGA-COAD","TCGA-UCS","TCGA-LUSC","TCGA-LGG","TCGA-HNSC","TCGA-TGCT","TARGET-CCSK","TCGA-THCA","TCGA-LIHC","TCGA-BLCA","TCGA-UCEC","TARGET-WT","TCGA-PRAD","TCGA-OV","TCGA-THYM","TCGA-CHOL","TCGA-STAD","TCGA-LUAD")

##### Store available tissue type codes and definitions in a list
Tissues <- as.list(setNames( c("Primary solid Tumor","Recurrent Solid Tumor","Primary Blood Derived Cancer - Peripheral Blood","Recurrent Blood Derived Cancer - Bone Marrow","Additional - New Primary","Metastatic","Additional Metastatic","Human Tumor Original Cells","Primary Blood Derived Cancer - Bone Marrow","Blood Derived Normal","Solid Tissue Normal","Buccal Cell Normal","EBV Immortalized Normal","Bone Marrow Normal","Control Analyte","Recurrent Blood Derived Cancer - Peripheral Blood","Cell Lines","Primary Xenograft Tissue","Cell Line Derived Xenograft Tissue"), c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,20,40,50,60,61)) )


##### Define the sample types
sampleType <- unlist(Tissues[ Tissue ])


##### Verify input parameters
if ( !is.na(args[4]) && !is.na(args[5]) && !is.na(args[5]) ) {

    ##### Verify project IDs
    if ( ProjectID %!in% ProjectIDs ) {

      cat( paste("Queried Project ID:", ProjectID, "is not available\n", sep=" ") )
      q()
    }

    ##### Verify tissue types
    if ( all(Tissue %!in% names(unlist(Tissues))) ) {

        cat(paste("\nSome of queried tissue types (", paste(Tissue, collapse = ", "), ") are not available!\n\n", sep=" "))
    }

} else {
    cat( "Some input parameters are mising or incorrect. Try again!" )
    q()
}


##### Set/create the project directory
if (file.exists(ProjectDir)){
	cat( paste("Folder with project [", ProjectName, "] already exists. The requested data will be stored in [", ProjectDir, "]\n", sep=" ") )
} else {
	dir.create(ProjectDir, recursive = TRUE);
}

##### Change working directory to the project workspace
setwd(ProjectDir)

##### Write used parameters into a file
write(args, file = "R_parameters.txt", append = FALSE, sep="\t")


#===============================================================================
#    Data query and download
#===============================================================================

##### Query DNA methylation data
query <- GDCquery( project = ProjectID, data.category = "DNA Methylation", platform = platform, sample.type = sampleType )


##### Download queried data. Use “client” method. Although the default "api" method is faster, but the data might get corrupted in the download, and it might need to be executed again. One can also set the "chunks.per.download" parameter. This will make the API method only download n files at a time. This may reduce the download problems when the data size is too large.
GDCdownload(query, method = "client" )


#===============================================================================
#    Data prepare
#===============================================================================

##### Prepare data matrix with the format described here https://gdc-docs.nci.nih.gov/Data/Bioinformatics_Pipelines/Methylation_LO_Pipeline/
data <- GDCprepare(query)


##### Use SummarizedExperiment package to extract the expression matrix from the SummarizedExperiment class used to store rectangular matrices of experimental results (sequencing and microarray experiments; http://bioconductor.org/packages/release/bioc/vignettes/SummarizedExperiment/inst/doc/SummarizedExperiment.html )

##### Get info about the features stored in the SummarizedExperiment DataFrame object
#rowData(data)
#rowRanges(data)

##### See header of the samples matrix
#colnames(SummarizedExperiment::colData(data))







##### NEED TO UPDATE FROM HERE

##### "shortLetterCode" is a column in the SummarizedExperiment::colData(data) matrix
TCGAvisualize_meanMethylation(data, groupCol = "shortLetterCode",filename = NULL)



##### Extract DNA methylation matrix
data_matrix <- assays(data)[[1]]

cat( paste( "Writing DNA methylation data [", paste(ProjectName, ".txt", sep=""), "] to [", ProjectDir, "]\n", sep=" ") )
write.table( prepare2write(data_matrix), file = paste(ProjectName, ".txt", sep="") ,sep="\t", row.names=FALSE )

##### Box plot of DNA methylation data per sample
pdf(paste(ProjectName, "_boxplot.pdf", sep=""), pointsize = 8 ,width = 0.1*ncol(data_matrix), height = 4)
par(mar=c(16,5,2,1))
boxplot( data_matrix,col="grey", las = 2, main = paste(ProjectName, " data", sep="") )
dev.off()



##### UPDATE UNTILL HERE





##### Extract samples information
samples_info <- colData(data)


cat( paste( "Writing samples information [", paste(ProjectName, "_samples.txt", sep=""), "] to [", ProjectDir, "]\n", sep=" ") )
write.table( samples_info, file = paste(ProjectName, "_samples.txt", sep="") ,sep="\t", row.names=FALSE )


#===============================================================================
#    Query and download associated clinical information
#===============================================================================

##### Query clinical data
query <- GDCquery( project = ProjectID, data.category = "Clinical", file.type = "xml" )

##### Download queried data. Use “client” method. Although the default "api" method is faster, but the data might get corrupted in the download, and it might need to be executed again. One can also set the "chunks.per.download" parameter. This will make the API method only download n files at a time. This may reduce the download problems when the data size is too large.
GDCdownload(query, method = "client")


#===============================================================================
#    Data prepare associated clinical information
#===============================================================================


clinical.patient <- GDCprepare_clinic(query, clinical.info = "patient")
rownames(clinical.patient) <- clinical.patient$bcr_patient_barcode


##### Create empty data frame for sample and clinical info
clinical.merged <- as.data.frame( setNames(replicate( sum(ncol(samples_info), ncol(clinical.patient) ), numeric(0), simplify = F), c(colnames(samples_info), colnames(clinical.patient))) )
clinical.present <- 0


##### Merge sample with clinical information
for (i in 1:nrow( samples_info ) ) {

    ##### Create empty data frame for clinical info
    clinical.2add <- as.data.frame( setNames(replicate(ncol(clinical.patient),numeric(0), simplify = F), colnames(clinical.patient)) )
    clinical.2add[1,] <- rep("", ncol(clinical.patient))

    ##### Scan the clincial information for each sample
    for (j in 1:nrow( clinical.patient ) ) {

        ##### Merge sample with clincial information if available for corresponding sample
        if ( samples_info$patient[i] == rownames(clinical.patient)[j] ) {

            clinical.2add <- clinical.patient[j,]
        }
    }
    clinical.merged <- rbind( clinical.merged, cbind( samples_info[i,], clinical.2add) )
}



cat( paste( "Writing samples and clinical information [", paste(ProjectName, "_clinical_info.txt", sep=""), "] to [", ProjectDir, "]\n", sep=" ") )
write.table( clinical.merged, file = paste(ProjectName, "_clinical_info.txt", sep="") ,sep="\t", row.names=FALSE )


##### Clear workspace
rm(list=ls())
##### Close any open graphics devices
graphics.off()
