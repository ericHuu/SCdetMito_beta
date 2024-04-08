#' SCQCone: Perform quality control (QC) for one sample/group single-cell RNA-seq data
#'
#' This function performs QC for a single-cell RNA-seq sample/group. It includes steps such as filtering cells based on library size,
#' mitochondrial content, and other relevant metrics. It generates various QC plots and tables to assess the quality of the data.
#'
#' @param seurat_obj A Seurat object containing single-cell RNA-seq data for a single sample
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
#' # devtools::install_github('satijalab/seurat-data')
#' # library(SeuratData)
#' # AvailableData()
#' # InstallData("pbmc3k")
#' # data("pbmc3k")
#' # pbmc3k$mitoRatio = PercentageFeatureSet(object = pbmc3k, pattern = "^MT-")/100
#' # seurat_obj <- SCQCone(pbmc3k)
#' # DO NOT RUN
#'
#' @seealso
#' Use this function to perform QC for a single-cell RNA-seq sample/group.
#'
#' @keywords internal
#'
#'

# Function for performing quality control on single-cell RNA-seq data
SCQCone <- function(
    # set seurat object
    seurat_obj,
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

    # If max_mito is set to "SCdetMito", calculate it using SCdetMito
    if (max_mito == "SCdetMito") {
        seurat_obj@meta.data$temped_sample <- rep(
            "A",
            ncol(seurat_obj)
        )
        max_mito <- SCdetMito(seurat_obj,
            by = "temped_sample",
            table_out = table_out,
            plot = plot
        )
        seurat_obj@meta.data <- subset(seurat_obj@meta.data,
            select = -temped_sample
        )
    } else if (max_mito != "SCdetMito") {
        max_mito <- max_mito
    }

    # Before QC
    # Generate summary table
    qc_summary <- data.frame(
        CellCounts = dim(seurat_obj)[2],
        GenesCounts_median = median(seurat_obj$nFeature_RNA),
        TranscriptCounts_median = median(seurat_obj$nCount_RNA),
        MitoRatio_median = median(seurat_obj$mitoRatio)
    )

    if (plot == TRUE) {
        seurat_obj@meta.data$temped_sample <- rep(
            "seurat_obj",
            ncol(seurat_obj)
        )
        SCQC_processedPlots(seurat_obj, by = "temped_sample", flag = "Checked")
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
        # %>%
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
        seurat_obj@meta.data$temped_sample <- rep(
            "seurat_obj",
            ncol(seurat_obj)
        )
        SCQC_processedPlots(seurat_obj, by = "temped_sample", flag = "Filtered")
    }

    # Print QC summary table
    message("QC Summary [after]:")
    print(qc_summary)

    # # Generate QC plots
    # if (plot == TRUE) {
    #     # Plot 1: Distribution of the number of genes per cell
    #     p1 <- ggplot2::ggplot(seurat_obj, ggplot2::aes(x = nFeature_RNA)) +
    #         ggplot2::geom_histogram(binwidth = 50, fill = "steelblue", color = "black") +
    #         ggplot2::theme_minimal() +
    #         ggplot2::labs(title = "Distribution of Genes per Cell", x = "Number of Genes", y = "Frequency")

    #     # Plot 2: Distribution of the total counts per cell
    #     p2 <- ggplot2::ggplot(seurat_obj, ggplot2::aes(x = nCount_RNA)) +
    #         ggplot2::geom_histogram(binwidth = 500, fill = "darkorange", color = "black") +
    #         ggplot2::theme_minimal() +
    #         ggplot2::labs(title = "Distribution of Counts per Cell", x = "Total Counts", y = "Frequency")

    #     # Plot 3: Distribution of mitochondrial gene ratio per cell
    #     p3 <- ggplot2::ggplot(seurat_obj, ggplot2::aes(x = percent.mt)) +
    #         ggplot2::geom_histogram(binwidth = 1, fill = "forestgreen", color = "black") +
    #         ggplot2::theme_minimal() +
    #         ggplot2::labs(title = "Distribution of Mitochondrial Ratio per Cell", x = "Percentage of Mitochondrial Genes", y = "Frequency")

    #     # Save the plots to a PDF file
    #     pdf("QC_Plots.pdf", width = 12, height = 5)
    #     print(p1)
    #     print(p2)
    #     print(p3)
    #     dev.off()
    # }

    # # Write out QC summary table
    # if (table_out == TRUE) {
    #     write.csv(qc_summary, "QC_Summary.csv", row.names = FALSE)
    # }

    return(seurat_obj)
}
