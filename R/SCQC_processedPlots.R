#' SCQC_processedPlots: Draw processed plots
#'
#' @param seurat_obj A Seurat object containing single-cell RNA-seq data
#' @param by Sample/group columns
#' @param flag Flag for out put files
#'
#' @return processed plots
#' @export
#'
#' @examples
#' # DO NOT RUN
#' # devtools::install_github('satijalab/seurat-data')
#' # library(SeuratData)
#' # AvailableData()
#' # InstallData("pbmc3k")
#' # data("pbmc3k")
#' # head(pbmc3k)
#' # library(Seurat)
#' # pbmc3k$mitoRatio <- PercentageFeatureSet(object = pbmc3k, pattern = "^MT-") / 100
#' # max(pbmc3k$mitoRatio)
#' # min(pbmc3k$mitoRatio)
#' # library(SCdetMito)
#' # dim(pbmc3k)
#' # pbmc3k$samples <- rep(
#' #    c("A", "B", "C"),
#' #    c(900, 900, 900)
#' #  )
#' # Idents(pbmc3k) <- pbmc3k$samples
#' # pbmc3k_A <- subset(pbmc3k,
#' #     idents = "A"
#' #  )
#' # qc： optional-1
#' #  qcpassed_pbmc3k_A <- SCQCone(pbmc3k_A,
#' #      max_mito = 0.06,
#' #      removeDouble = FALSE
#' #  )
#' #  dim(qcpassed_pbmc3k_A)
#' # qc： optional-2
#' #  qcpassed_pbmc3k_A <- SCQCone(pbmc3k_A,
#' #      max_mito = 0.06,
#' #      removeDouble = FALSE
#' #  )
#' #  dim(qcpassed_pbmc3k_A)
#' # qc： optional-3
#' #  qcpassed_pbmc3k_A <- SCQCone(pbmc3k_A,
#' #      max_mito = "SCdetMito",
#' #      removeDouble = FALSE
#' #  )
#' #  dim(qcpassed_pbmc3k_A)
#'
#' @seealso
#' Use this function to generate processed plots
#'
#' @keywords internal
#'
#'
SCQC_processedPlots <- function(seurat_obj,by,flag) {

    # head(seurat_obj@meta.data)
    metadata <- seurat_obj@meta.data
    colnames(metadata)
    metadata$cells <- rownames(metadata)

    ########################
    # library(RColorBrewer)
    qual_col_pals <- RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == "qual", ]
    # 处理后有73种差异还比较明显的颜色，基本够用
    col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    set.seed(6)
    cols <- sample(col_vector, length(unique(metadata[[`by`]])))
    cols

    # 1 对 cell counts 绘图，正常的样本数在6000-12000左右
    # metadata %>%
        ggplot2::ggplot(metadata,ggplot2::aes(x = .data[[`by`]], fill = .data[[`by`]])) +
        ggplot2::scale_fill_manual(values = cols) +
        ggplot2::geom_bar(width = 0.6) +
        ggplot2::theme_bw() +
        ggplot2::theme(plot.margin = ggplot2::unit(rep(1.5, 4), "cm")) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1)) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")) +
        ggplot2::ggtitle("Cells number of each sample") +
        Seurat::NoLegend() +
        ggplot2::labs(x = "Samples", y = "Counts", size = 2)
    ggplot2::ggsave(paste0(flag,"_Ncells_in_each_sample.pdf"), width = 8, height = 5)
    while (!is.null(dev.list())) dev.off()

    # 2 UMI counts (transcripts) per cell
    # metadata %>%
        ggplot2::ggplot(metadata,ggplot2::aes(x = nCount_RNA, fill = .data[[`by`]])) +
        ggplot2::scale_fill_manual(values = cols) +
        ggplot2::geom_density(alpha = 0.67) +
        ggplot2::theme_bw() +
        ggplot2::theme(plot.margin = ggplot2::unit(rep(1.5, 4), "cm")) +
        ggplot2::scale_x_log10() +
        ggplot2::ggtitle("Cell density per UMI number") +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")) +
        ggplot2::guides(fill = ggplot2::guide_legend(title = "Samples")) +
        ggplot2::labs(x = "UMI numbers", y = "Cell density", size = 2) 
        # ggplot2::geom_vline(xintercept = 500, col = "tomato", size = 0.6)
        ggplot2::ggsave(paste0(flag,"_Cell_density_per_UMI.pdf"), width = 8, height = 5)
    while (!is.null(dev.list())) dev.off()

    # 3 Genes detected per cell
    # metadata %>%
        ggplot2::ggplot(metadata,ggplot2::aes(x = nFeature_RNA, fill = .data[[`by`]])) +
        ggplot2::scale_fill_manual(values = cols) +
        ggplot2::geom_density(alpha = 0.67) +
        ggplot2::theme_bw() +
        ggplot2::theme(plot.margin = ggplot2::unit(rep(1.5, 4), "cm")) +
        ggplot2::scale_x_log10() +
        ggplot2::ggtitle("Cell density per gene number") +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")) +
        ggplot2::guides(fill = ggplot2::guide_legend(title = "Samples")) +
        ggplot2::labs(x = "Gene numbers", y = "Cell density", size = 2)
    ggplot2::ggsave(paste0(flag,"_Cell_density_per_Gene.pdf"), width = 8, height = 5)
    while (!is.null(dev.list())) dev.off()

    # metadata %>%
        ggplot2::ggplot(metadata,ggplot2::aes(x = .data[[`by`]], y = log10(nFeature_RNA), fill = .data[[`by`]])) +
        ggplot2::scale_fill_manual(values = cols) +
        ggplot2::geom_boxplot(outlier.size = 0.5) +
        ggplot2::guides(fill = "none") +
        ggplot2::theme_bw() +
        ggplot2::theme(plot.margin = ggplot2::unit(rep(1.5, 4), "cm")) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1)) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")) +
        ggplot2::ggtitle("Genes per cell in each sample") +
        ggplot2::labs(x = "Samples", y = "log10(Gene number)", size = 2)
    ggplot2::ggsave(paste0(flag,"_gene_per_cell.pdf"), width = 8, height = 5)
    while (!is.null(dev.list())) dev.off()

    # 3.4 UMI v.s. genes detected
    # metadata %>%
        ggplot2::ggplot(metadata,ggplot2::aes(x = nCount_RNA, y = nFeature_RNA, color = mitoRatio)) +
        ggplot2::geom_point(size = 0.5) +
        ggplot2::scale_color_gradient(low = "lightgray", high = "tomato") +
        ggplot2::stat_smooth(method = lm) +
        ggplot2::scale_x_log10() +
        ggplot2::scale_y_log10() +
        ggplot2::theme_bw() +
        ggplot2::theme(plot.margin = ggplot2::unit(rep(1.5, 4), "cm")) +
        ggplot2::facet_wrap(~ .data[[`by`]])
    ggplot2::ggsave(paste0(flag,"_gene_vs_umi-mitoRatio.pdf"), width = 12, height = 10)
    while (!is.null(dev.list())) dev.off()

    # 5 mitochondrial counts ratio
    # metadata %>%
        ggplot2::ggplot(metadata,ggplot2::aes(x = mitoRatio, fill =  .data[[`by`]])) +
        ggplot2::scale_fill_manual(values = cols) +
        ggplot2::geom_density(alpha = 0.5, size = 0.1) +
        ggplot2::guides(fill = ggplot2::guide_legend(title = "Samples")) +
        # scale_x_log10()+
        ggplot2::theme_bw() +
        ggplot2::theme(plot.margin = ggplot2::unit(rep(1.5, 4), "cm")) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")) +
        ggplot2::ggtitle("Cells distribution by MitoRatio") +
        ggplot2::labs(x = "MitoRatio", y = "Density", size = 2) 
        # ggplot2::geom_vline(xintercept = 0.2) # 根据数据来判断
        ggplot2::ggsave(paste0(flag,"_mitochondrial.pdf"), width = 9, height = 5)
    while (!is.null(dev.list())) dev.off()

    # 3.6 Complexity
    metadata$log10GenesPerUMI <- log10(seurat_obj$nFeature_RNA) / log10(seurat_obj$nCount_RNA)
    # metadata %>%
        ggplot2::ggplot(metadata,ggplot2::aes(x = log10GenesPerUMI, fill = .data[[`by`]])) +
        ggplot2::scale_fill_manual(values = cols) +
        ggplot2::theme_bw() +
        ggplot2::theme(plot.margin = ggplot2::unit(rep(1.5, 4), "cm")) +
        ggplot2::geom_density(alpha = 0.5, size = 0.2) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")) +
        ggplot2::ggtitle("Genes complexity") +
        ggplot2::labs(x = "log10(Genes Per UMI)", y = "Density", size = 2) +
        # ggplot2::geom_vline(xintercept = 0.8) +
        ggplot2::guides(fill = ggplot2::guide_legend(title = "Samples"))
    ggplot2::ggsave(paste0(flag,"_genes_Complexity.pdf"), width = 9, height = 5)
    while (!is.null(dev.list())) dev.off()

    # # make a place for rRNA (to be added)
    # # 3.7 rRNA
    # metadata %>%
    #     ggplot2::ggplot(ggplot2::aes(x = percent.r_genes, fill = .data[[`by`]])) +
    #     ggplot2::scale_fill_manual(values = cols) +
    #     ggplot2::theme_bw() +
    #     ggplot2::theme(plot.margin = ggplot2::unit(rep(1.5, 4), "cm")) +
    #     ggplot2::geom_density(alpha = 0.5, size = 0.2) +
    #     ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")) +
    #     ggplot2::ggtitle("  ") +
    #     ggplot2::labs(x = "percentage of ribosome RNA", y = "Density", size = 2) +
    #     # geom_vline(xintercept = 0.8)+
    #     ggplot2::guides(fill = ggplot2::guide_legend(title = "Samples"))
    # ggplot2::ggsave(paste0(flag,"_percent.r_genes.pdf"), width = 9, height = 5)
    # while (!is.null(dev.list())) dev.off()

    # metadata %>%
    #     ggplot2::ggplot(ggplot2::aes(x = nCount_RNA, y = nFeature_RNA, color = percent.r_genes)) +
    #     ggplot2::geom_point(size = 0.5) +
    #     ggplot2::scale_color_gradient(low = "lightgray", high = "tomato") +
    #     ggplot2::stat_smooth(method = lm) +
    #     ggplot2::scale_x_log10() +
    #     ggplot2::scale_y_log10() +
    #     ggplot2::theme_bw() +
    #     ggplot2::theme(plot.margin = unit(rep(1.5, 4), "cm")) +
    #     ggplot2::geom_vline(xintercept = 500) +
    #     ggplot2::geom_hline(yintercept = 250) +
    #     ggplot2::facet_wrap(~.data[[`by`]])
    # ggplot2::ggsave(paste0(flag,"_gene_vs_umi-percent.r_genes.pdf"), width = 12, height = 10)
    # while (!is.null(dev.list())) dev.off()
}
