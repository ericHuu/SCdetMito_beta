#' check_seu: checking if mitoRatio (or any colnumn) is present in your Seurat object
#'
#' @param seurat_obj Provide your Seurat object
#' @param check Provide a column name [default: 'mitoRatio']
#'
#' @return Checked Seurat object
#' @export
#'
#' @examples
#'# DO NOT RUN 
#'# devtools::install_github('satijalab/seurat-data')
#'# library(SeuratData)
#'# AvailableData()
#'# InstallData("pbmc3k")
#'# data("pbmc3k")
#'# pbmc3k$mitoRatio = PercentageFeatureSet(object = pbmc3k,pattern = "^MT-")/100
#'# check_seu(pbmc3k,"mitoRatio")
#'# DO NOT RUN
#'
#'
#' @seealso
#' Use this function before performing mito-cutoff change point detection with SCdetMito.
#' 
#' @keywords internal
#'
#' 
check_seu <- function(seurat_obj, 
check="mitoRatio") {
      # Check if the specified column is present in the meta.data of the Seurat object,
      # especially 'mitoRatio'
    if (!check %in% colnames(seurat_obj@meta.data)) {
        stop("The specified column '", 
        check, 
        "' is not present in the Seurat object. \n 
        Please prepare the data before checking!!!")
    }
      # If the specified column is 'mitoRatio' and contains values greater than 1, 
      # and assume percentage values and divide by 100
    if (check=="mitoRatio"){
        if (max(seurat_obj[[check]])>1){
            seurat_obj[[check]] = seurat_obj[[check]]/100
        }
    }
    return(seurat_obj)
}
