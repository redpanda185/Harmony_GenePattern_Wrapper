library(harmony)
library(dplyr)
library(Seurat)
library(optparse)
library(cowplot)
library(xfun)

parser = OptionParser()
parser <- add_option(parser, c("--file_list"), 
                     help = "list of files to read", 
                     default = "NO_FILE")
parser <- add_option(parser, c("--output_file"), type = 'character',
                     help = "Name of the file to save the data to", 
                     default = "harmonized_data")

print('================================================')
args <- parse_args(parser)
print('Parameters used:')
print(args)
print('================================================')

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
    row.names(readTable) <- readTable$symbol
    readTable[1] <- NULL
    
    # print(readTable[1:6,1:6])
    
    panels[[i]] = as(as.matrix(readTable), "dgCMatrix")
    names(panels)[[i]] = as.character(name)
    i = i + 1
  }
}
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
}
save_it <- function(object, fileName){
  saveRDS(object, file = fileName)
  print("Saved file!")
  return(object)
}
output <- run_harmony(panels)
save_it(output, paste(args$output_file, '.rds', sep = ''))
