#01 SCdetMito_beta
Implements a comprehensive set of functions for detecting change points in the mitochondrial ratio (mitoRatio) across single-cell samples. SCdetMito is specifically designed for the analysis of mitochondrial content variations in single-cell RNA sequencing data. It offers methods to identify significant change points and visualize patterns of mitoRatio alterations, making it particularly useful for exploring mitochondrial dynamics in single-cell studies. The package contributes to data quality control in single-cell RNA-seq experiments.

#02 install by
install_github("ericHuu/SCdetMito_beta")

#03 Load
library(SCdetMito_beta)

#04 QC for One sample 
#04-1 InstallData("pbmc3k")
data("pbmc3k")
library(Seurat)
pbmc3k$mitoRatio <- PercentageFeatureSet(object = pbmc3k, pattern = "^MT-") / 100
# generate test samples
dim(pbmc3k)
pbmc3k$samples <- rep(
    c("A", "B", "C"),
    c(900, 900, 900)
)

## 挑选出样本 A 并进行 QC
Idents(pbmc3k) <- pbmc3k$samples
pbmc3k_A <- subset(pbmc3k,
    idents = "A"
)

# qc： optional-1
# 自定义max_mito 6%, 对双细胞过滤
qcpassed_pbmc3k_A <- SCQCone(pbmc3k_A,
    max_mito = 0.06,
    removeDouble = FALSE
)
dim(qcpassed_pbmc3k_A) # 查看过滤后数据结构

# qc： optional-2
# 自定义max_mito 6%, 不对双细胞过滤
qcpassed_pbmc3k_A <- SCQCone(pbmc3k_A,
    max_mito = 0.06,
    removeDouble = FALSE
)
dim(qcpassed_pbmc3k_A) # 查看过滤后数据结构

# qc： optional-3
# 使用 SCdetMito 计算最大 mitoRatio，并移除双重峰值
qcpassed_pbmc3k_A <- SCQCone(pbmc3k_A,
    max_mito = "SCdetMito",
    removeDouble = FALSE
)
dim(qcpassed_pbmc3k_A) # 查看过滤后数据结构
