library(devtools)
library(harmony)
library(dplyr)
library(Seurat)
library(optparse)
library(cowplot)
library(purrr)
library(xfun)

parser = OptionParser()
parser <- add_option(parser, c("--file_list"), 
                     help = "Matrix with coordinates for each cell (row) along many PCs (columns).", 
                     default = "NO_FILE")

print('================================================')
args <- parse_args(parser)
print('Parameters used:')
print(args)
print('================================================')

# RUN_PREPROCESS <<- FALSE
# 
# if(tolower(args.data_preprocess) == "yes"){
#   RUN_PREPROCESS <<- TRUE
# }

print(args$file_list)
con <- file(args$file_list, open = "r")
lines = readLines(con)
panels = list(NA)
i = 1
for (line in lines) {
  print("About to read")
  print(line)
  if (file_ext(line) == "rds") {
    RDSinput = readRDS(line)
    name = tail(strsplit(line, "/")[[1]], 1)
    name <- gsub("\\.rds$", "", name)
    print(paste0("Using ", name, " as the name of the column"))
    panels[[i]] = RDSinput
    names(panels)[[i]] = as.character(name)
    i = i + 1
  } else {
    name = tail(strsplit(line, "/")[[1]], 1)
    print(paste0("Using ", name, " as the name of the column"))
    readTable = read.table(line, sep = "\t", header = TRUE, )
    
    # readTable = readTable[1:500,1:501] write.table(readTable, file =
    # paste('small_500x500',name,sep='_'), sep = '\t',row.names = F, quote=F)
    
    row.names(readTable) <- readTable$symbol
    readTable[1] <- NULL
    
    # print(readTable[1:6,1:6])
    
    panels[[i]] = as(as.matrix(readTable), "dgCMatrix")
    names(panels)[[i]] = as.character(name)
    i = i + 1
  }
}

# do_preprocess(args){
#   data_mat1 <- args$data_mat1
#   rownames(data_mat1) <- data_mat1[, 1]
#   colnames(data_mat1) <- data_mat1[1, ]
#   data_mat1 <- data_mat1[2:nrow(data_mat1), 2:ncol(data_mat1)]
#   data_mat2 <- args$data_mat2
#   rownames(data_mat2) <- data_mat2[, 1]
#   colnames(data_mat2) <- data_mat2[1, ]
#   data_mat2 <- data_mat2[2:nrow(data_mat2), 2:ncol(data_mat2)]
#   cell_lines <- cbind(data_mat1, data_mat2)
#   obj <- CreateSeuratObject(counts = cell_lines, project = "harmony", min.cells = 3, min.genes = 200) %>%
#     Seurat::NormalizeData(verbose = FALSE) %>%
#     FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
#     ScaleData(verbose = FALSE) %>% 
#     RunPCA(pc.genes = pbmc@var.genes, npcs = 20, verbose = FALSE)
#   obj$celltype <- c(rep("Type 1", ncol(data_mat1)), rep("Type 2", ncol(data_mat2)))  #Change nammes of cell types to something more appropriate or something that can be input
#   return(obj)
# }
run_harmony <- function(datalist){
  data <- datalist[[1]]
  metadata_list <- NULL
  for(x in 1:length(datalist)){
    if(x != 1){
      data <- merge(data, datalist[[x]], merge.dr = "pca", merge.data = TRUE)
    }
    metadata_list <- c(metadata_list, rep(names(datalist)[x], ncol(datalist[[x]])))
  }
  data$celltype <- metadata_list
  data <- ScaleData(data)
  harmonizedData <- data %>%
    RunHarmony(group.by = "celltype", plot_convergence = TRUE)
  # p1 <- DimPlot(object = harmonizedData, reduction = "harmony", pt.size = .1, group.by = "celltype")
  # p2 <- VlnPlot(object = harmonizedData, features = "harmony_1", group.by = "celltype", pt.size = .1)
  # plot_grid(p1,p2)
  return(harmonizedData)
  # if(RUN_PREPROCESS){
  #   seuratObj <- do_prepreocess(args)
  #   harmonizedobj <- seuratObj %>%
  #     RunHarmony(group.by = "celltype", plot_convergence = TRUE)
  #   options(repr.plot.height = 5, repr.plot.width = 12)
  #   p1 <- DimPlot(object = harmonizedobj, reduction = "harmony", pt.size = .1, group.by = "celltype")
  #   p2 <- VlnPlot(object = harmonizedobj, features = "harmony_1", group.by = "celltype", pt.size = .1)
  #   plot_grid(p1,p2)
  # }
  # else{
  #   data_mat <- read.table(args$data_mat)
  #   meta_dat <- read.table(args$meta_data)
  #   harmonized_pcs <- harmony::HarmonyMatrix(
  #     data_matrix  = data_mat,       # Matrix with coordinates for each cell (row) along many PCs (columns)
  #     meta_data = meta_data, # Dataframe with information for each cell (row)
  #     vars_use  = args$vars_use, # Column in meta_data that defines dataset for each cell
  #     do_pca    = args$do_pca      # Since we are providing PCs, do not run PCA
  #   )
  # }
}
run_harmony(panels)
