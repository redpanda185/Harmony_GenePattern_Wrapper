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
parser <- add_option(parser, c("--cell_types"),
                     help = "Names of the cell types present in the dataset, separated by commas",
                     default = NULL)
parser <- add_option(parser, c("--colors"),
                     help = "The colors you want to use to label your cell types, separated by commas",
                     default = NULL)
# =================================================
# Advanced options
parser <- add_option(parser, c("--reduction"),
                     help = "Name of dimension reduction to use. Default is PCA",
                     default = 'pca')
parser <- add_option(parser, c("--dims.use"),
                     help = "Which PCA dimensions to use for Harmony. By default, use all",
                     default = NULL)
parser <- add_option(parser, c("--theta"),
                     help = "Diversity clustering penalty parameter. Specify for each variable 
                     in group.by.vars. Default theta=2. theta=0 does not encourage any 
                     diversity. Larger values of theta result in more diverse clusters",
                     default = NULL)
parser <- add_option(parser, c("--lambda"),
                     help = "Ridge regression penalty parameter. Specify for each variable
                     in group.by.vars. Default lambda=1. Lambda must be strictly positive.
                     Smaller values result in more aggressive correction.",
                     default = NULL)
parser <- add_option(parser, c("--sigma"),
                     help = "Width of soft kmeans clusters. Default sigma=0.1. Sigma scales
                     the distance from a cell to cluster centroids. Larger values of sigma result
                     in cells assigned to more clusters. Smaller values of sigma make soft
                     kmeans cluster approach hard clustering.",
                     default = 0.1)
parser <- add_option(parser, c("--nclust"),
                     help = "Number of clusters in model. nclust=1 equivalent to simple
                     linear regression.",
                     default = NULL)
parser <- add_option(parser, c("--tau"),
                     help = "Protection against overclustering small datasets with large ones.
                     tau is the expected number of cells per cluster.",
                     default = 0)
parser <- add_option(parser, c("--block.size"),
                     help = "What proportion of cells to update during clustering.
                     Between 0 to 1, default 0.05. Larger values may be faster but less accurate",
                     default = 0.05)
parser <- add_option(parser, c("--max.iter.harmony"),
                     help = "Maximum number of iterations that Harmony will run",
                     type = 'integer',
                     default = 10)
parser <- add_option(parser, c("--max.iter.cluster"),
                     help = "Maximum number of rounds to run clustering at each
                     round of Harmony.",
                     type = 'integer',
                     default = 20)
parser <- add_option(parser, c("--epsilon.cluster"),
                     help = "Convergence tolerance for clustering round of Harmony.
                     Set to -Inf to never stop early.",
                     default = 1e-5)
parser <- add_option(parser, c("--epsilon.harmony"),
                     help = "Convergence tolerance for Harmony. Set to -Inf to
                     never stop early.",
                     default = 1e-4)
parser <- add_option(parser, c("--plot_convergence"),
                     help = "Whether to print the convergence plot of the
                     clustering objective function. TRUE to plot, FALSE to suppress. This can be
                     useful for debugging.",
                     default = FALSE)
parser <- add_option(parser, c("--verbose"),
                     help = "Whether to print progress messages. TRUE to print, FALSE to
                     suppress.",
                     default = TRUE)
parser <- add_option(parser, c("--reference_values"),
                     help = "Defines reference dataset(s). Cells
                     that have batch variables values matching reference_values will not be moved",
                     default = NULL)
parser <- add_option(parser, c("--reduction.save"),
                     help = "Keyword to save Harmony reduction. Useful if you want
                     to try Harmony with multiple parameters and save them as e.g.
                     'harmony_theta0', 'harmony_theta1', 'harmony_theta2'",
                     default = "harmony")
parser <- add_option(parser, c("--assay.use"),
                     help = "Which assay to run PCA on if no PCA present?",
                     default = NULL)
parser <- add_option(parser, c("--project.dim"),
                     help = "Project dimension reduction loadings. Default TRUE.",
                     default = TRUE)


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
#Checks if the number of cell types matches the number of files given.
if(!is.null(namelist)){
  namelist <- strsplit(namelist, split = ',')
  if(length(namelist[[1]]) != length(lines)) stop("Number of cell types does not match number of files")
  lapply(namelist, trimws)
}

colorlist = args$colors
#Checks if the number of colors matches the number of files given.
if(!is.null(colorlist)){
  colorlist <- strsplit(colorlist, split = ',')
  if(length(colorlist[[1]]) != length(lines)) stop("Number of colors does not match number of cell types")
}

#Reads in the files from file_list
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

#Runs Harmony on the given data and outputs the processed data in a Seurat object.
run_harmony <- function(datalist, args){
  data <- datalist[[1]]
  metadata_list <- NULL
  #Merges all the Seurat objects into one object, while keeping the pca data. 
  for(x in 1:length(datalist)){
    if(x != 1){
      data <- merge(data, datalist[[x]], merge.dr = "pca", merge.data = TRUE)
    }
    metadata_list <- c(metadata_list, rep(names(datalist)[x], ncol(datalist[[x]])))
  }
  #Sets the metadata of the Seurat object to describe the cell type from which the data originates.
  data$celltype <- metadata_list
  #Runs all of the pre-processing steps prior to running Harmony
  data <- data %>%
    Seurat::NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(pc.genes = pbmc@var.genes, npcs = 20, verbose = FALSE)
  harmonizedData <- data %>%
    RunHarmony(group.by = "celltype", 
               reduction = args$reduction,
               dims.use = args$dims.use,
               theta = args$theta,
               lambda = args$lambda,
               sigma = args$sigma,
               nclust = args$nclust,
               tau = args$tau,
               block.size = args$block.size,
               max.iter.harmony = args$max.iter.harmony,
               max.iter.cluster = args$max.iter.cluster,
               epsilon.cluster = args$epsilon.cluster,
               epsilon.harmony = args$epsilon.harmony,
               plot_convergence = args$plot_convergence,
               verbose = args$verbose,
               reference_values = args$reference_values,
               reduction.save = args$reduction.save,
               assay.use = args$assay.use,
               project.dim = args$project.dim
               )
  return(harmonizedData)
}

#Saves a Seurat object to an RDS file
save_it <- function(object, fileName){
  saveRDS(object, file = fileName)
  print("Saved file!")
  return(object)
}

#Produces a scatterplot of the harmonized data and saves it as a png file
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

output <- run_harmony(panels, args)
save_it(output, paste(args$output_file, '.rds', sep = ''))
visualize(output, colorlist[[1]])
