#' SCdetMito: Do mitoRatio change point detection
#'
#' @param seurat_obj A Seurat object containing single-cell RNA-seq data
#' @param mitoRatio The column name for mitoRatio in the Seurat object
#' @param by The target variable for neighbor comparison would be the 'samples' column, which would be used to conduct counts distribution and tests for each sample. The result would be a mitoRation cutoff for those sample. The target variable colnumn could be 'group' or 'treatment' either.
#' @param bin_width 0.01(default), slide window width between every two cutoff, any value less than 0.1 could be used
#' @param min_cut 0.01(default), minimum cutoff during the detection process, any value less than 0.1 could be used
#' @param max_cut 1(default), maximum cutoff during the detection process, any value less than 1 could be used
#' @param table_out Set to TRUE or FALSE(default) for table generation
#' @param plot Set to TRUE or FALSE(default) for plots generation
#' @param diff_num Cell changed-counts, which at least filtered by a slide cut off
#' @param ... Additional parameters
#'
#' @return mitoRatio cutoff value and (optional) some processed tables and plots
#' @export
#'
#' @examples
#' # DO NOT RUN
#' # devtools::install_github('satijalab/seurat-data')
#' # library(SeuratData)
#' # AvailableData()
#' # InstallData("pbmc3k")
#' # data("pbmc3k")
#' # pbmc3k
#' # library(Seurat)
#' # pbmc3k$mitoRatio = PercentageFeatureSet(object = pbmc3k,
#' #                                         pattern = "^MT-")/100
#' # head(pbmc3k)
#' #
#' # library(SCdetMito)
#' #
#' # take "orig.ident" as one sample
#' # mitoRatio_cut_off = SCdetMito(pbmc3k,
#' #                                by="orig.ident",
#' #                                table_out = T)
#' #
#' #
#' # generate 3 samples for test
#' # dim(pbmc3k)
#' # pbmc3k$samples = rep(c("A","B","C"),
#' #                      c(900,900,900))
#' # mitoRatio_cut_off = SCdetMito(pbmc3k,
#' #                               by="samples",
#' #                               table_out = T)
#' # DO NOT RUN
#'
#' 
#' @seealso
#' Use this function to identify mitoRatio change points in your single-cell RNA-seq data.
#'
#' @keywords internal
#'
#' 
SCdetMito <- function(
    seurat_obj,
    mitoRatio = "mitoRatio",
    by = "sample",
    bin_width = 0.01,
    min_cut = 0.01,
    max_cut = 1,
    table_out = FALSE,
    plot = TRUE,
    diff_num = 5,
    ...) {
    # Print messages using cat() for clearer console output
    message("Processing the mito-cutoff change point detecting...")

    # Check and print progress
    message("Checking the infos of the target and mitoRatio...")
    seurat_obj <- check_seu(seurat_obj, by)
    seurat_obj <- check_seu(seurat_obj, mitoRatio)
    message("Checking the infos of the target and mitoRatio OK!")

    raw_seurat_obj <- seurat_obj
    ## test, and debugs 
    # print(by)
    # print(raw_seurat_obj[[by]])
    sample_counts_by_mito <- table(raw_seurat_obj[[by]])
    sample_counts_by_mito <- as.data.frame(sample_counts_by_mito)
    # print(colnames(sample_counts_by_mito))
    colnames(sample_counts_by_mito) <- c(`by`, "raw_counts")
    # head(sample_counts_by_mito)

    ## Perform cutoff detection
    # 1 0.99 0.98 ...
    j <- max_cut / bin_width
    k <- min_cut / bin_width

    # Subset the raw Seurat object and calculate sample counts by mitoRatio
    for (i in seq(j, k)) {
        # Subset the raw Seurat object based on nUMI and mitoRatio
        raw_seurat_obj_subset <- subset(
            x = raw_seurat_obj,
            subset = mitoRatio < bin_width * i
        )

        # Calculate sample counts for the subset
        sample_counts <- as.data.frame(table(raw_seurat_obj_subset[[by]]))

        # Add the sample counts to a new data frame
        sample_counts_by_mito <- cbind(sample_counts_by_mito, sample_counts$Freq)

        # Rename the new columns with more accurate names
        colnames(sample_counts_by_mito)[ncol(sample_counts_by_mito)] <- paste0(i * bin_width)
    }
    ## take a look at counts distrubution by mitoRatio
    print(head(sample_counts_by_mito))

    # Write out counts distribution table
    if (table_out == TRUE) {
        write.csv(sample_counts_by_mito,
         "processed_temp_cell_num.csv", 
         quote = FALSE, 
         row.names = FALSE
         )
    }


    # Reshape the sample_counts_by_mito data frame into a long format for plotting
    data <- reshape2::melt(sample_counts_by_mito)
    ## take a look at reshaped-counts before plots-drawing
    print(head(data))

    # Create a ggplot object with the counts data, with each sample represented by a different color
    if (plot == TRUE) {
        p <- ggplot2::ggplot(
            data,
            ggplot2::aes(x = variable, 
            y = value, 
            group = .data[[`by`]], 
            color = .data[[`by`]])
        ) +
            ggplot2::geom_line(linetype = "dashed", 
            size = 0.8) +
            ggplot2::geom_point(size = 1) +
            ggplot2::theme(
                panel.background = ggplot2::element_blank(),
                axis.line = ggplot2::element_line(colour = "black"),
                panel.border = ggplot2::element_rect(colour = "black", 
                fill = NA),
                axis.text.x = ggplot2::element_text(angle = 45, 
                hjust = 1, 
                size = 6)
            ) +
            ggplot2::xlab(`by`) +
            ggplot2::ylab("Counts") +
            ggplot2::scale_color_manual(values = rainbow(length(unique(data[[by]]))))

        # Save the plot to a PDF file
        ggplot2::ggsave("CountsDistributionBy_Mito.pdf", 
        width = 12, 
        height = 5)

        # Check device and turn off
        while (!is.null(dev.list())) dev.off()
    }


    # Create an empty data frame to store the results and an empty list to store the plots
    results <- data.frame()
    plots <- list()

    # Loop through each row (sample) in the sample_counts_by_mito data frame
    for (i in 1:nrow(sample_counts_by_mito)) {
        # Extract the counts for the current sample and remove the first column (which contains the sample name)
        n <- sample_counts_by_mito[i, ][, 2:ncol(sample_counts_by_mito)]

        # Create an empty vector to store the positions of the identified mito-cutoffs
        cut <- c()

        # Loop through each position in the counts data frame (excluding the first and last positions)
        for (k in seq(2, (ncol(n) - 1))) {
            # Check if the difference between the current position and the adjacent positions is greater than 5
            if ((n[k - 1] - n[k] > diff_num) || (n[k] - n[k + 1] > diff_num)) {
                # If the difference is greater than 5, create a data frame with the counts for the current position and the adjacent positions
                a <- data.frame(count_bef = c(n[k - 1], n[k]), 
                count_af = c(n[k - 1] - n[k], n[k] - n[k + 1]))

                # Perform a chi-squared test on the data frame to test for a significant difference between the counts
                p <- chisq.test(a, 
                correct = F)

                # Adjust the p-value using the false discovery rate (FDR) method
                p_fdr <- p.adjust(p$p.value, 
                method = "fdr")

                # Print the position and p-value to the console
                print(paste0(k, "    ", p$p.value))

                # Add the position to the cut vector
                cut <- c(cut, k)

                # Add the results for the current sample and position to the results data frame
                results <- rbind(results, 
                data.frame(by = (sample_counts_by_mito[[`by`]])[i], 
                cutoff = as.numeric(colnames(n)[k]), 
                position = k,
                p_value = p$p.value, 
                p_fdr = p_fdr))
            }
        }

        # Reshape the counts data frame into a long format for plotting
        data <- reshape2::melt(n)

        # Create a ggplot object with the counts data and add vertical lines at the identified mito-cutoff positions
        p <- ggplot2::ggplot(
            data,
            ggplot2::aes(x = variable, 
            y = value, 
            group = 1)
        ) +
            ggplot2::geom_line(ggplot2::aes(color = "Data"), 
            size = 1, 
            linetype = "dashed") +
            ggplot2::geom_point(ggplot2::aes(color = "Data"), 
            size = 1) +
            ggplot2::geom_vline(xintercept = cut, 
            linetype = 2, 
            linewidth = 0.68, 
            color = "tomato") +
            ggplot2::theme(
                panel.background = ggplot2::element_blank(),
                axis.line = ggplot2::element_line(colour = "black"),
                panel.border = ggplot2::element_rect(colour = "black", 
                fill = NA),
                axis.text.x = ggplot2::element_text(angle = 45, 
                hjust = 1)
            ) +
            ggplot2::xlab("Cut-offs") +
            ggplot2::ylab("Counts") +
            ggplot2::scale_color_manual(values = c("Data" = "steelblue"))

        # Add the sample name to the plot title and add a subtitle
        p <- p + ggplot2::ggtitle((sample_counts_by_mito[[`by`]])[i]) +
            ggplot2::labs(subtitle = "Mito-change-point-detect")

        # Add the plot to the plots list
        plots[[i]] <- p
    }


    # Combine all plots into a single plot grid using the grid.arrange() function
    p_final <- gridExtra::grid.arrange(grobs = plots, 
    ncol = 1)
    i <- i * 0.5
    # Save the final plot and the results data frame to files
    if (plot == TRUE) {
        ggplot2::ggsave(paste0("your-mito-change-point-detect", ".pdf"),
            p_final,
            width = 12, 
            height = 5 * i, 
            limitsize = F
        )
        # Check device and turn off
        while (!is.null(dev.list())) dev.off()
    }


    # Create an empty list to store significant points for each sample
    significant_points <- list()

    # Loop through each unique sample in the results data frame
    for (i in unique(results$by)) {
        print(i)
        # Filter the results data frame to include only the current sample and significant points with FDR-adjusted p-value <= 0.05
        sig_points <- results[which(results$by == i & results$p_fdr <= 0.05), ]
        sig_points <- sig_points$cutoff
        # Add the significant points for the current sample to the significant_points list
        significant_points[[i]] <- sig_points
    }

    # Combine all significant points into a single vector
    all_points <- unlist(significant_points)

    # Count the frequency of each significant point
    freq <- table(all_points)

    # Identify significant points that occur in at least 50% of samples
    significant_freq <- freq[freq >= length(unique(results[[`by`]])) * .5]
    significant_freq <- significant_freq[order(names(significant_freq))]

    # Print the maximum significant point as the final mito-cutoff
    message("Chosed mito-cutoff: ", max(as.numeric(names(significant_freq))), "\n")
    my_final_cutoff <- max(as.numeric(names(significant_freq)))

    return(my_final_cutoff)

    colnames(results)[1] <- `by`
    # write out significant mitoRatios to table
    if (table_out == TRUE) {
        write.csv(results, 
        "mito_change_point_results.csv", 
        row.names = FALSE)
    }
}
