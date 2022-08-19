library(devtools)
library(harmony)
library(dplyr)
library(Seurat)
library(optparse)
library(cowplot)

parser = OptionParser()
parser <- add_option(parser, c("--data_preprocess"),
                     help = "Does the data need to be preprocessed",
                     deafault = "YES")
parser <- add_option(parser, c("--data_mat1"), 
                     help = "Matrix with coordinates for each cell (row) along many PCs (columns).", 
                     default = "NO_FILE")
parser <- add_option(parser, c("--data_mat2"), 
                      help = "Matrix with coordinates for each cell (row) along many PCs (columns).", 
                      default = "NO_FILE")
parser <- add_option(parser, c("--do_pca"), 
                     type = 'boolean',
                     help = "Since we are providing PCs, do not run PCA.", 
                     default = "NO_FILE")

print('================================================')
args <- parse_args(parser)
print('Parameters used:')
print(args)
print('================================================')

RUN_PREPROCESS <<- FALSE

if(tolower(args.data_preprocess) == "yes"){
  RUN_PREPROCESS <<- TRUE
}

do_preprocess(args){
  data_mat1 <- args$data_mat1
  rownames(data_mat1) <- data_mat1[, 1]
  colnames(data_mat1) <- data_mat1[1, ]
  data_mat1 <- data_mat1[2:nrow(data_mat1), 2:ncol(data_mat1)]
  data_mat2 <- args$data_mat2
  rownames(data_mat2) <- data_mat2[, 1]
  colnames(data_mat2) <- data_mat2[1, ]
  data_mat2 <- data_mat2[2:nrow(data_mat2), 2:ncol(data_mat2)]
  cell_lines <- cbind(data_mat1, data_mat2)
  obj <- CreateSeuratObject(counts = cell_lines, project = "harmony", min.cells = 3, min.genes = 200) %>%
    Seurat::NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
    ScaleData(verbose = FALSE) %>% 
    RunPCA(pc.genes = pbmc@var.genes, npcs = 20, verbose = FALSE)
  obj$celltype <- c(rep("Type 1", ncol(data_mat1)), rep("Type 2", ncol(data_mat2)))  #Change nammes of cell types to something more appropriate or something that can be input
  return(obj)
}
run_harmony <- function(args){
  if(RUN_PREPROCESS){
    seuratObj <- do_prepreocess(args)
    harmonizedobj <- seuratObj %>%
      RunHarmony(group.by = "celltype", plot_convergence = TRUE)
    options(repr.plot.height = 5, repr.plot.width = 12)
    p1 <- DimPlot(object = harmonizedobj, reduction = "harmony", pt.size = .1, group.by = "celltype")
    p2 <- VlnPlot(object = harmonizedobj, features = "harmony_1", group.by = "celltype", pt.size = .1)
    plot_grid(p1,p2)
  }
  else{
    data_mat <- read.table(args$data_mat)
    meta_dat <- read.table(args$meta_data)
    harmonized_pcs <- harmony::HarmonyMatrix(
      data_matrix  = data_mat,       # Matrix with coordinates for each cell (row) along many PCs (columns)
      meta_data = meta_data, # Dataframe with information for each cell (row)
      vars_use  = args$vars_use, # Column in meta_data that defines dataset for each cell
      do_pca    = args$do_pca      # Since we are providing PCs, do not run PCA
    )
  }
}