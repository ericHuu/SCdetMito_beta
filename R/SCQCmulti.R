#' SCQCmulti: Perform quality control (QC) for one sample/group single-cell RNA-seq data
#'
#' This function performs QC for a single-cell RNA-seq sample/group. It includes steps such as filtering cells based on library size,
#' mitochondrial content, and other relevant metrics. It generates various QC plots and tables to assess the quality of the data.
#' @param seurat_obj A Seurat object containing single-cell RNA-seq data for a single sample
#' @param by The target variable for neighbor comparison would be the 'samples' column, which would be used to conduct counts distribution and tests for each sample. The result would be a mitoRation cutoff for those sample. The target variable colnumn could be 'group' or 'treatment' either.
#' @param mode How to QC, "split" means QC for each sample/group and then merge; "all" means QC as a sample/group
#' @param min_genes Minimum number of genes expressed to consider a cell [default: 200]
#' @param max_genes Maximum number of genes expressed to consider a cell [default: Inf]
#' @param min_counts Minimum total counts per cell to consider a cell [default: 500]
#' @param max_counts Maximum total counts per cell to consider a cell [default: Inf]
#' @param min_mito Minimum percentage of mitochondrial genes per cell to consider a cell [default: 0]
#' @param max_mito Maximum percentage of mitochondrial genes per cell to consider a cell [default: 20]
#' @param removeDouble Do doublets cell or not [default: TRUE]
#' @param plot Set to TRUE to generate QC plots [default: TRUE]
#' @param table_out Set to TRUE to generate QC summary table [default: TRUE]
#' @param ... Additional parameters (not used)
#'
#' @return QC summarys, QC-passed seurat object, and, if specified, QC plots
#'
#' @examples
#' # DO NOT RUN
#' # # devtools::install_github('satijalab/seurat-data')
#' # library(SeuratData)
#' # AvailableData()
#' # # InstallData("pbmc3k")
#' # data("pbmc3k")
#' # head(pbmc3k)
#' # library(Seurat)
#' # pbmc3k$mitoRatio <- PercentageFeatureSet(object = pbmc3k, pattern = "^MT-") / 100
#' # max(pbmc3k$mitoRatio)
#' # min(pbmc3k$mitoRatio)
#' # # 加载 SCdetMito 包
#' # library(SCdetMito)
#' # dim(pbmc3k)
#' # pbmc3k$samples <- rep(
#' #     c("A", "B", "C"),
#' #     c(900, 900, 900)
#' # )
#' # head(pbmc3k)
#' # # qc： optional-1
#' # qcpassed_pbmc3k <- SCQCmulti(pbmc3k,
#' #     by = "samples",
#' #     mode = "all",
#' #     max_mito = "SCdetMito",
#' #     removeDouble = FALSE
#' # )
#' # dim(qcpassed_pbmc3k)
#' # # qc： optional-2
#' # qcpassed_pbmc3k <- SCQCmulti(pbmc3k,
#' #     by = "samples",
#' #     mode = "all",
#' #     max_mito = "SCdetMito",
#' #     removeDouble = T
#' # )
#' # dim(qcpassed_pbmc3k)
#' # # qc： optional-3
#' # qcpassed_pbmc3k <- SCQCmulti(pbmc3k,
#' #     by = "samples",
#' #     mode = "split",
#' #     max_mito = "SCdetMito",
#' #     removeDouble = FALSE
#' # )
#' # DO NOT RUN
#'
#' @seealso
#' Use this function to perform QC for a single-cell RNA-seq sample/group.
#'
#' @keywords internal
#'
#'

# Function for performing quality control on single-cell RNA-seq data
SCQCmulti <- function(
    # set seurat object
    seurat_obj,
    by,
    mode = c("split", "all"),
    # set used colnumns
    nFeature_RNA = "nFeature_RNA", # nGene
    nCount_RNA = "nCount_RNA", # nUMI
    mitoRatio = "mitoRatio", # percent.mt
    # set for filters
    min_genes = 200,
    max_genes = Inf,
    min_counts = 500,
    max_counts = Inf,
    min_mito = 0,
    max_mito = "SCdetMito", # or set a float range from 0-1
    removeDouble = TRUE,
    plot = TRUE,
    table_out = TRUE,
    ...) {
    message("Performing quality control for single-cell RNA-seq data... ")

    by <- by
    # print(by)
    mode <- match.arg(mode)

    # Check if specified columns exist in Seurat object
    seurat_obj <- check_seu(
        seurat_obj,
        nFeature_RNA
    )
    seurat_obj <- check_seu(
        seurat_obj,
        nCount_RNA
    )
    seurat_obj <- check_seu(
        seurat_obj,
        mitoRatio
    )

    if (mode == "all") {
        # print(mode)
        # If max_mito is set to "SCdetMito", calculate it using SCdetMito
        if (max_mito == "SCdetMito") {
            max_mito <- SCdetMito(seurat_obj,
                by = by,
                table_out = table_out,
                plot = plot
            )
        } else if (max_mito != "SCdetMito") {
            max_mito <- max_mito
        }
        # print(max_mito)

        # Before QC
        # Generate summary table
        qc_summary <- data.frame(
            CellCounts = dim(seurat_obj)[2],
            GenesCounts_median = median(seurat_obj$nFeature_RNA),
            TranscriptCounts_median = median(seurat_obj$nCount_RNA),
            MitoRatio_median = median(seurat_obj$mitoRatio)
        )

        if (plot == TRUE) {
            meta.data <- seurat_obj@meta.data
            meta.data["temped_sample"] <- meta.data[`by`]
            seurat_obj@meta.data <- meta.data
            SCQC_processedPlots(seurat_obj, by = "temped_sample", flag = "M-Checked")
        }

        # Print summary table
        message("QC Summary [before]:")
        print(qc_summary)

        # Display main QC parameters
        message("mitoRatio cutoff was set to: ", max_mito)
        message("Other QC indicators: ")
        message("min_genes: ", min_genes)
        message("max_genes: ", max_genes)
        message("min_counts: ", min_counts)
        message("max_counts: ", max_counts)

        # Filtering cells based on library size, mitochondrial content, and other relevant metrics
        message("Perform cell filtering ...")
        seurat_obj_F <- subset(
            x = seurat_obj,
            subset = (nFeature_RNA >= min_genes) &
                (nFeature_RNA <= max_genes) &
                (nCount_RNA >= min_counts) &
                (nCount_RNA <= max_counts) &
                (mitoRatio >= min_mito) &
                (mitoRatio <= max_mito)
        )

        message("Perform double-cell filtering ...")
        # If removeDouble is TRUE, identify and remove doublets
        if (removeDouble == TRUE) {
            seurat_obj <- seurat_obj_F
            seurat_obj <- NormalizeData(seurat_obj) # default "LogNormalize"
            seurat_obj <- FindVariableFeatures(
                object = seurat_obj,
                selection.method = "vst",
                nfeatures = 2000,
                verbose = F
            )
            seurat_obj <- ScaleData(object = seurat_obj)
            seurat_obj <- RunPCA(
                object = seurat_obj,
                features = VariableFeatures(object = seurat_obj)
            )
            seurat_obj <- RunUMAP(seurat_obj,
                dims = 1:30
            )
            pc.num <- 1:30
            sweep.res.list <- DoubletFinder::paramSweep_v3(seurat_obj,
                PCs = pc.num,
                sct = F
            )
            sweep.stats <- DoubletFinder::summarizeSweep(sweep.res.list,
                GT = FALSE
            )
            bcmvn <- DoubletFinder::find.pK(sweep.stats)
            pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)]
            #  %>%
            pK_bcmvn <- as.character(pK_bcmvn)
            pK_bcmvn <- as.numeric(pK_bcmvn)
            DoubletRate <- ncol(seurat_obj) / 5000 * 0.039
            message(paste0("PLS cking DoubletRate: ", DoubletRate))
            nExp_poi <- round(DoubletRate * ncol(seurat_obj))

            ## Identification of doublets
            seurat_obj <- DoubletFinder::doubletFinder_v3(seurat_obj,
                PCs = pc.num,
                pN = 0.25,
                pK = pK_bcmvn,
                nExp = nExp_poi,
                reuse.pANN = F,
                sct = F
            )
            colnames(seurat_obj@meta.data)[ncol(seurat_obj@meta.data)] <- "double_info"
            # Remove singlets
            dou_rm <- subset(seurat_obj,
                subset = double_info == "Singlet"
            )
            seurat_obj <- dou_rm
        } else if (removeDouble == FALSE) {
            seurat_obj <- seurat_obj_F
        }
        # After QC
        # Generate QC summary table
        qc_summary <- data.frame(
            CellCounts = dim(seurat_obj)[2],
            GenesCounts_median = median(seurat_obj$nFeature_RNA),
            TranscriptCounts_median = median(seurat_obj$nCount_RNA),
            MitoRatio_median = median(seurat_obj$mitoRatio)
        )

        if (plot == TRUE) {
            meta.data <- seurat_obj@meta.data
            meta.data["temped_sample"] <- meta.data[`by`]
            seurat_obj@meta.data <- meta.data
            SCQC_processedPlots(seurat_obj, 
            by = "temped_sample", 
            flag = "M-Filtered")
        }

        # Print QC summary table
        message("QC Summary [after]:")
        # print(qc_summary)

        return(seurat_obj)
    } else if (mode == "split") {
        seurat_obj_tmp <- seurat_obj
        if (plot == TRUE) {
            meta.data <- seurat_obj_tmp@meta.data
            meta.data["temped_sample"] <- meta.data[`by`]
            seurat_obj_tmp@meta.data <- meta.data
            SCQC_processedPlots(seurat_obj_tmp, 
            by = "temped_sample", 
            flag = "M-Checked")
        }
        ns <- unique(seurat_obj_tmp@meta.data[[`by`]])
        nn <- length(unique(seurat_obj_tmp@meta.data[[`by`]]))
        y <- c()
        for (i in seq(nn)) {
            Idents(seurat_obj_tmp) <- seurat_obj_tmp[[`by`]]
            seurat_obj_tmp_selcted <- subset(seurat_obj_tmp,
                idents = ns[i]
            )
            qcpassed_seurat_obj_tmp_selcted <- SCQCone(seurat_obj_tmp_selcted,
                nFeature_RNA = nFeature_RNA,
                nCount_RNA = nCount_RNA,
                mitoRatio = mitoRatio,
                min_genes = 200,
                max_genes = Inf,
                min_counts = 500,
                max_counts = Inf,
                min_mito = 0,
                max_mito = "SCdetMito",
                removeDouble = TRUE,
                plot = F,
                table_out = F,
            )

            assign(paste0("a", i), 
            qcpassed_seurat_obj_tmp_selcted)
            if (i != 1) {
                y <- c(y, c(assign(paste0("a", i), 
                qcpassed_seurat_obj_tmp_selcted)))
            }
        }
        qcpassed_seurat_obj <- merge(
            x = a1,
            y = y
        )
        return(qcpassed_seurat_obj)
        if (plot == TRUE) {
            meta.data <- qcpassed_seurat_obj@meta.data
            meta.data["temped_sample"] <- meta.data[`by`]
            qcpassed_seurat_obj@meta.data <- meta.data
            SCQC_processedPlots(qcpassed_seurat_obj, 
            by = "temped_sample", 
            flag = "M-Filtered")
        }
    }
}
