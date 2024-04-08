# SCdetMito_beta: An R Package for Detecting Optimal mitoRatio in Single-Cell RNA-Seq Data Quality Control
Version: 1.0
Authors@R: person("Silu Hu", email = "erichu121@foxmail.com", role = c("aut", "cre"))

# 01 Introduction
Implements a comprehensive set of functions for detecting change points in the mitochondrial ratio (mitoRatio) across single-cell samples. SCdetMito is specifically designed for the analysis of mitochondrial content variations in single-cell RNA sequencing data. It offers methods to identify significant change points and visualize patterns of mitoRatio alterations, making it particularly useful for exploring mitochondrial dynamics in single-cell studies. The package contributes to data quality control in single-cell RNA-seq experiments.

# 02 install by
install_github("ericHuu/SCdetMito_beta")

# 03 Load
library(SCdetMito_beta)

# 04 QC for One sample 
## 04-1 take pbmc3k as example
InstallData("pbmc3k")
data("pbmc3k")
library(Seurat)
pbmc3k$mitoRatio <- PercentageFeatureSet(object = pbmc3k, pattern = "^MT-") / 100
## 04-2 Generate test samples
dim(pbmc3k)
pbmc3k$samples <- rep(
    c("A", "B", "C"),
    c(900, 900, 900)
)
## 04-3 Select A for QC
Idents(pbmc3k) <- pbmc3k$samples
pbmc3k_A <- subset(pbmc3k,
    idents = "A"
)
### qc： optional-1， set max_mito to 6, removeDouble = TRUE
qcpassed_pbmc3k_A <- SCQCone(pbmc3k_A,
    max_mito = 0.06,
    removeDouble = TRUE
)
### Check data after QC
dim(qcpassed_pbmc3k_A)
### qc： optional-2, set max_mito to 6%, removeDouble = FALSE
qcpassed_pbmc3k_A <- SCQCone(pbmc3k_A,
    max_mito = 0.06,
    removeDouble = FALSE
)
### Check data after QC
dim(qcpassed_pbmc3k_A)
### qc： optional-3, by SCdetMito (7%), removeDouble = TRUE
qcpassed_pbmc3k_A <- SCQCone(pbmc3k_A,
    max_mito = "SCdetMito",
    removeDouble = TRUE
)
![image](https://github.com/ericHuu/SCdetMito_beta/blob/main/img/your-mito-change-point-detect.png)
### Check data after QC
dim(qcpassed_pbmc3k_A) 

# 05 QC for multiple samples
## 05-1 take pbmc3k as example
InstallData("pbmc3k")
data("pbmc3k")
library(Seurat)
pbmc3k$mitoRatio <- PercentageFeatureSet(object = pbmc3k, pattern = "^MT-") / 100
## 05-2 Generate test samples
dim(pbmc3k)
pbmc3k$samples <- rep(
    c("A", "B", "C"),
    c(900, 900, 900)
)
head(pbmc3k)
## 05-3 QC for 3 samples
### qc： optional-1, set max_mito as 6%, removeDouble = TRUE (or FALSE)
qcpassed_pbmc3k <- SCQCmulti(pbmc3k,
    by = "samples",
    mode = "all",
    max_mito = "SCdetMito",
    removeDouble = TRUE
)
### Check data after QC
dim(qcpassed_pbmc3k)

### qc： optional-2, by SCdetMito, removeDouble = TRUE, mode = "all"
qcpassed_pbmc3k <- SCQCmulti(pbmc3k,
    by = "samples",
    mode = "all",
    max_mito = "SCdetMito",
    removeDouble = T
)
### Check data after QC
dim(qcpassed_pbmc3k)

### qc： optional-3, by SCdetMito, removeDouble = FALSE, mode = "split"
qcpassed_pbmc3k <- SCQCmulti(pbmc3k,
    by = "samples",
    mode = "split",
    max_mito = "SCdetMito",
    removeDouble = FALSE
)
### Check data after QC
dim(qcpassed_pbmc3k)


