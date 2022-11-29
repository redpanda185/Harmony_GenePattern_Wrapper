library(harmony)
library(dplyr)
library(Seurat)
library(optparse)
library(cowplot)
library(xfun)

parser = OptionParser()
parser <- add_option(parser, c("--file_list"), 
                     help = "List of Seurat files to read", 
                     default = "NO_FILE")
parser <- add_option(parser, c("--output_file"), type = 'character',
                     help = "Name of the file to save the data to", 
                     default = "harmonized_data")
parser <- add_option(parser, c("--max_iterations"),
                     help = "Maximum number of iterations that Harmony will run")
parser <- add_option(parser, c("--cell_types"),
                     help = "Names of the cell types present in the dataset, separated by commas",
                     default = NULL)
parser <- add_option(parser, c("--colors"),
                     help = "The colors you want to use to label your cell types, separated by commas",
                     default = NULL)
# parser <- add_option(parser, c("--number of cells"),
#                      help = "Number of cells per cell types (must be in the same order)")

print('================================================')
args <- parse_args(parser)
print('Parameters used:')
print(args)
print('================================================')

print(args$file_list)
con <- file(args$file_list, open = "r")
lines = readLines(con)
panels = list(NA)

namelist = args$cell_types
if(!is.null(namelist)){
  namelist <- strsplit(namelist, split = ',')
  if(length(namelist[[1]]) != length(lines)) stop("Number of cell types does not match number of files")
  lapply(namelist, trimws)
}
colorlist = args$colors
if(!is.null(colorlist)){
  colorlist <- strsplit(colorlist, split = ',')
  if(length(colorlist[[1]]) != length(lines)) stop("Number of colors does not match number of cell types")
}

max_iter = args$max_iterations
if(typeof(max_iter) == "character"){
  max_iter <- as.numeric(max_iter)
}
i = 1

for (line in lines) {
  print("About to read")
  print(line)
  if (file_ext(line) == "rds") {
    RDSinput = readRDS(line)
    name = NULL
    if(is.null(args$cell_types)){
      name = tail(strsplit(line, "/")[[1]], 1)
      name <- gsub("\\.rds$", "", name)
    }
    else{
      name = namelist[[1]][i]
    }
    print(paste0("Using ", name, " as the name of the column"))
    panels[[i]] = RDSinput
    names(panels)[[i]] = as.character(name)
    i = i + 1
  } else {
    name = tail(strsplit(line, "/")[[1]], 1)
    print(paste0("Using ", name, " as the name of the column"))
    readTable = read.table(line, sep = "\t", header = TRUE, )
    row.names(readTable) <- readTable$symbol
    readTable[1] <- NULL
    
    # print(readTable[1:6,1:6])
    
    panels[[i]] = as(as.matrix(readTable), "dgCMatrix")
    names(panels)[[i]] = as.character(name)
    i = i + 1
  }
}
run_harmony <- function(datalist, max_iterations){
  data <- datalist[[1]]
  metadata_list <- NULL
  for(x in 1:length(datalist)){
    if(x != 1){
      data <- merge(data, datalist[[x]], merge.dr = "pca", merge.data = TRUE)
    }
    metadata_list <- c(metadata_list, rep(names(datalist)[x], ncol(datalist[[x]])))
  }
  data$celltype <- metadata_list
  data <- data %>%
    Seurat::NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(pc.genes = pbmc@var.genes, npcs = 20, verbose = FALSE)
  harmonizedData <- data %>%
    RunHarmony(group.by = "celltype", plot_convergence = TRUE, max.iter.harmony = max_iterations)
  return(harmonizedData)
}
save_it <- function(object, fileName){
  saveRDS(object, file = fileName)
  print("Saved file!")
  return(object)
}
visualize <- function(harmonyObj, colors){
  print("Making plot!")
  if(is.null(colors)){
    print("Using default colors to make plot!")
    p1 <- DimPlot(object = harmonyObj, reduction = "harmony", pt.size = .1, group.by = "celltype")
  }
  else{
    print("Plot colors:")
    lapply(colors, print)
    p1 <- DimPlot(object = harmonyObj, reduction = "harmony", pt.size = .1, group.by = "celltype", cols = colors)
  }
  png(filename = paste(args$output_file, '_plot' ,'.png', sep = ''), width = 2000, height = 1600, res = 400)
  print(paste0("Using ", paste(args$output_file, '_plot' ,'.png', sep = ''), " as file name for the plot."))
  plot(p1)
  dev.off()
  print("Plot made and saved!")
}
output <- run_harmony(panels, max_iter)
save_it(output, paste(args$output_file, '.rds', sep = ''))
visualize(output, colorlist[[1]])
