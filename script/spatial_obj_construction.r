library(Seurat)
library(SeuratData)
library(dplyr)
library("argparse")
p <- ArgumentParser(description="construction")
p$add_argument("--path", type="character",required=T)
p$add_argument("-o", "--outdir", type="character",required=T)
p$add_argument("-p", "--prefix", type="character", default="gene")
p$add_argument("--png", type="character",required=T)
opt <- p$parse_args()
if( !file.exists(opt$outdir) ){
  if( !dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE) ){
    stop(paste("failed to create the path: ",opt$outdir,sep=""))
  }
}
expr <- Seurat::Read10X(opt$path, cell.column = 1)
obj <- Seurat::CreateSeuratObject(counts = expr,
                              assay = 'Spatial',
                               min.cells=1,
                            min.features=1)
barcode_pos_path <- paste0(opt$path,'/barcodes_pos.tsv.gz')
barcode_pos <- read.table(gzfile(barcode_pos_path),header = F) %>%
    dplyr::rename(Barcode = V1 , pos_w = V2, pos_h = V3)
barcode_pos <- barcode_pos %>% dplyr::filter(., Barcode %in% rownames(obj@meta.data))
#make spatial coord file for seurat S4 class
coord <- data.frame(tissue = 1,
                      row = barcode_pos$pos_h,
                      col = barcode_pos$pos_w,
                      imagerow = barcode_pos$pos_h,
                      imagecol = barcode_pos$pos_w)
rownames(coord) <- barcode_pos$Barcode
png <- png::readPNG(opt$png)
zoom_scale <- 1.02316090204888
sample1 <-  new(Class = "VisiumV1",
                  image = png,
                  scale.factors = Seurat::scalefactors(zoom_scale, 100, zoom_scale, zoom_scale),
                  coordinates = coord,
                  spot.radius = 0.0045,
                  assay = 'Spatial',
                  key = "sample1_")
obj@images <- list(sample1 = sample1)
DefaultAssay(obj) <- "Spatial"
message(DefaultAssay(obj))
saveRDS(obj, file = paste(opt$outdir,"/",opt$prefix,".rds", sep=""))



