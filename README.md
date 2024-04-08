# 01 SCdetMito_beta
Implements a comprehensive set of functions for detecting change points in the mitochondrial ratio (mitoRatio) across single-cell samples. SCdetMito is specifically designed for the analysis of mitochondrial content variations in single-cell RNA sequencing data. It offers methods to identify significant change points and visualize patterns of mitoRatio alterations, making it particularly useful for exploring mitochondrial dynamics in single-cell studies. The package contributes to data quality control in single-cell RNA-seq experiments.

# 02 install by
install_github("ericHuu/SCdetMito_beta")

# 03 Load
library(SCdetMito_beta)

# 04 QC for One sample 
## 04-1 InstallData("pbmc3k")
data("pbmc3k")
library(Seurat)
pbmc3k$mitoRatio <- PercentageFeatureSet(object = pbmc3k, pattern = "^MT-") / 100
## 04-2 generate test samples
dim(pbmc3k)
pbmc3k$samples <- rep(
    c("A", "B", "C"),
    c(900, 900, 900)
)
## 04-3 select A for QC
Idents(pbmc3k) <- pbmc3k$samples
pbmc3k_A <- subset(pbmc3k,
    idents = "A"
)
### qc： optional-1， set max_mito to 6, removeDouble = TRUE
qcpassed_pbmc3k_A <- SCQCone(pbmc3k_A,
    max_mito = 0.06,
    removeDouble = TRUE
)
### check data after QC
dim(qcpassed_pbmc3k_A)
### qc： optional-2, set max_mito to 6%, removeDouble = FALSE
qcpassed_pbmc3k_A <- SCQCone(pbmc3k_A,
    max_mito = 0.06,
    removeDouble = FALSE
)
### check data after QC
dim(qcpassed_pbmc3k_A)
### qc： optional-3, by SCdetMito, removeDouble = TRUE
qcpassed_pbmc3k_A <- SCQCone(pbmc3k_A,
    max_mito = "SCdetMito",
    removeDouble = TRUE
)
### check data after QC
dim(qcpassed_pbmc3k_A) 
