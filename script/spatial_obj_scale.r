library(Seurat)
library(SeuratData)
library(dplyr)
library("argparse")
p <- ArgumentParser(description="heatmap")
p$add_argument("-r", "--rds", type="character",required=T)
p$add_argument("-o", "--outdir", type="character",required=T)
p$add_argument("-p", "--prefix", type="character", default="gene")

opt <- p$parse_args()
if( !file.exists(opt$outdir) ){
  if( !dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE) ){
    stop(paste("dir.create failed: outdir=",opt$outdir,sep=""))
  }
}

obj <- readRDS(opt$rds)
DefaultAssay(obj) <- "Spatial"
message(DefaultAssay(obj))
obj <- NormalizeData(obj, normalization.method = "LogNormalize",scale.factor = 1e6)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
obj <- ScaleData(obj, features = rownames(obj), verbose = TRUE, vars.to.regress = NULL)

saveRDS(obj, file = paste(opt$outdir,"/",opt$prefix,"_scale.rds", sep=""))






