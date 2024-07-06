library(Seurat)
library("argparse")

p <- ArgumentParser(description="edit meta data")

p$add_argument("-r", "--rds", type="character",required=T)
p$add_argument("-o", "--outdir", type="character",required=T)
p$add_argument("-p", "--prefix", type="character", required=T)
p$add_argument("--albst", type="character",required=T)
p$add_argument("--start", type="double",default=0)
p$add_argument("--end", type="double",default=2)

opt <- p$parse_args()
if( !file.exists(opt$outdir) ){
  if( !dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE) ){
    stop(paste(" failed to create path: ",opt$outdir,sep=""))
  }
}

obj <- NULL
obj <- readRDS(opt$rds)
meta <- NULL
meta <- obj@meta.data
head(meta)


meta$ID <- row.names(meta)


albst <- read.table(opt$albst, header = T, sep = " ")

meta <- merge(meta, albst, by = "ID", all.x = T, sort = F)

row.names(meta) <- meta$ID

start <- NULL
end <- NULL
step <- NULL
thresholds <- NULL



start <- opt$start
end <- opt$end
step <- 0.05
thresholds <- seq(start, end, by = step)

for (i in 1:(length(thresholds) - 1)) {
	if (i == 1){
	meta$pos[meta$ALBST >= thresholds[i] & meta$ALBST <= thresholds[i + 1]] <- thresholds[i]
	}else{
	meta$pos[meta$ALBST > thresholds[i] & meta$ALBST <= thresholds[i + 1]] <- thresholds[i]
    }
}

unique(meta$pos)
obj@meta.data <- meta
head(obj@meta.data)
saveRDS(obj, file = paste0(opt$outdir,"/",opt$prefix,".rds"))








