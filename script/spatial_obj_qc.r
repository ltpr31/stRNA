library("argparse")
library(Seurat)
library(future)
options(future.globals.maxSize = 50 * 1024^3)

parser <- ArgumentParser(description='Spatial spots quality control')
parser$add_argument( "--rds", type="character",required=F)
parser$add_argument( "--nGene.max", type="integer", default=NULL)
parser$add_argument( "--nGene.min", type="integer", default=NULL)
parser$add_argument( "--nUMI.min", type="integer", default=NULL)
parser$add_argument( "--nUMI.max", type="integer", default=NULL)
parser$add_argument( "--nspot.min", type="integer", default=5)
parser$add_argument( "-o", "--outdir", type="character", default=getwd())
parser$add_argument( "-p", "--prefix", type="character", required=T)

opt <- parser$parse_args()

if( !file.exists(opt$outdir) ){
  if( !dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE) ){
    stop(paste("failed to create the path:",opt$outdir,sep=""))
  }
}

obj<- readRDS(opt$rds)

if (is.null(opt$nUMI.max)) {
	nUMI.max <- Inf
}else{
	nUMI.max=opt$nUMI.max
}

if (is.null(opt$nUMI.min)) {
	nUMI.min <- -Inf
}else{
	nUMI.min=opt$nUMI.min
}

if (is.null(opt$nGene.max)) {
	nGene.max <- Inf
}else{
	nGene.max=opt$nGene.max
}
if (is.null(opt$nGene.min)) {
	nGene.min <- -Inf
}else{
	nGene.min=opt$nGene.min
}

cat("before filter: ", dim(obj), "\n")
obj <- subset(obj, subset = nFeature_Spatial > nGene.min & nFeature_Spatial < nGene.max  & nCount_Spatial < nUMI.max & nCount_Spatial > nUMI.min )

meta <- obj@meta.data

value_counts <- table(meta$pos)
value_counts_df <- as.data.frame(table(meta$pos))
filt <- value_counts_df[value_counts_df$Freq >= opt$nspot.min,]
obj <- subset(obj, subset = pos %in% filt$Var1)

cat("after filter: ", dim(obj), "\n")
saveRDS(obj, file = paste0(opt$outdir,"/",opt$prefix,".rds"))


