

library(Seurat)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(pheatmap)
library(tidyverse)
library(ggplotify)
library("argparse")
p <- ArgumentParser(description="heatmap")

p$add_argument("-l", "--list", type="character", required=T)
p$add_argument("-r", "--rds", type="character",required=T)
p$add_argument("-o", "--outdir", type="character",required=T)
p$add_argument("-p", "--prefix", type="character", default="gene")
p$add_argument("-H", "--height", type="double", default=0.2)
p$add_argument("-W", "--width", type="double", default=0.2)
p$add_argument("--cluster", action="store_true", default=FALSE)
p$add_argument("--heatmap", type="double", default=30)
p$add_argument("--num", type="double", default=20)
p$add_argument("--miny", type="double", default=0)
p$add_argument("--fit", type="double", default=60)
p$add_argument("--gnum", type="double", default=60)
p$add_argument("--holdall", action="store_true", default=FALSE)
p$add_argument("--start", type="double", default=0)
p$add_argument("--end", type="double", default=2)

opt <- p$parse_args()
if( !file.exists(opt$outdir) ){
  if( !dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE) ){
    stop(paste("dir.create failed: outdir=",opt$outdir,sep=""))
  }
}

obj <- readRDS(opt$rds)
DefaultAssay(obj) <- "Spatial"
message(DefaultAssay(obj))
list <- read.csv(opt$list, sep = "\t", header = TRUE)
list$ID <- gsub("_","-",list$ID)
IDs <- list$ID
IDs <- unique(IDs)
gene_cell_exp <- AverageExpression(obj,
                                   features = IDs,
                                   group.by = "pos",
                                   slot = 'data')
merge <- as.data.frame(gene_cell_exp$Spatial)
write.table(merge,paste(opt$outdir,"/",opt$prefix,"_spatial_dataframe.txt", sep=""), sep = "\t", quote = F, col.names = T, row.names = T)
if(opt$holdall){
  
  merge$ID <- row.names(merge)
  merge$ID <- gsub("-","_",merge$ID)
  IDs <- list$ID
  IDs <- unique(IDs)
  IDs <- as.data.frame(IDs)
  
  colnames(IDs)[colnames(IDs) == 'IDs'] <- 'ID'
  IDs$ID <- gsub("-","_",IDs$ID)
  
  merge <- merge(IDs, merge, by = "ID", all.x = T, sort = F)
  merge[is.na(merge)] <- 0
  row.names(merge) <- merge$ID
  merge$ID <- NULL
}else{
  merge <- merge[rowSums(merge != 0) > 0, ]
}

nw <- dim(merge)[2]
w = (opt$width*nw)
nh <- dim(merge)[1]
h = (opt$height*nh) + 8
merge[is.na(merge)] <- 0
if(opt$cluster & dim(merge)[1]>=2){
  
  p1 <- pheatmap(merge, 
                 color = colorRampPalette(c( "white", "#862B0D"))(100), 
                 scale = "none", 
                 cluster_col = F,
                 cluster_row = T,
                 cellwidth = 7,
                 cellheight = 1,
                 angle_col = 90,
                 border = 1,
                 border_color = "black",
                 legend_border_color = "black",
                 show_rownames = FALSE,
                 fontsize_row = 2,
                 fontsize_col = 5,
                 treeheight_row=10,
                 annotation_legend_colors = "black",
                 main = "spatial_heatmap",
                 breaks = seq(0, opt$heatmap,length.out = 101))
  
  
}else{
  p1 <- pheatmap(merge, 
                 color = colorRampPalette(c( "white", "#862B0D"))(100), 
                 scale = "none", 
                 cluster_col = F,
                 cluster_row = F,
                 cellwidth = 7,
                 cellheight = 2,
                 angle_col = 90,
                 border = 1,
                 border_color = "black",
                 legend_border_color = "black",
                 show_rownames = FALSE,
                 fontsize_row = 2,
                 fontsize_col = 5,
                 treeheight_row=10,
                 annotation_legend_colors = "black",
                 main = "spatial_heatmap",
                 breaks = seq(0, opt$heatmap,length.out = 101))
  
}
p1_grob <- as.grob(p1)


gene_cell_exp <- AverageExpression(obj,
                                   group.by = "pos",
                                   slot = 'data')
all_gene <- as.data.frame(gene_cell_exp$Spatial)
all_gene[is.na(all_gene)] <- 0
all_gene <- all_gene[rowSums(all_gene != 0) > 0, ]
all_gene <- t(all_gene)
all_gene <- as.data.frame(all_gene)
all_gene$all_num <- apply(all_gene, 1, function(x) sum(x > opt$num))
all_gene$ALBST <- row.names(all_gene)
all_gene <- all_gene[,c("ALBST","all_num")]


data <- merge
data <- t(data)
data <- as.data.frame(data)
data$num <- apply(data, 1, function(x) sum(x > opt$num))
data$ALBST <- row.names(data)
data <- data[,c("ALBST","num")]



data <- merge(data, all_gene, all = T, by = "ALBST", sort = F)

data$num <- as.numeric(data$num)
data$all_num <- as.numeric(data$all_num)

data$rate <- 1000*data$num/data$all_num

data$num <- data$rate

df_bar <- data
df_bar$ALBST <- as.numeric(df_bar$ALBST)
class(df_bar)
# calculation of fit value
fit <- mgcv::gam(num ~ s(ALBST, bs = "cs"), data = df_bar)
df_bar$fit <- predict(fit, newdata = df_bar)
df_bar <- arrange(df_bar, ALBST)
df_bar$fit2 <-  ifelse(df_bar$fit > opt$fit, opt$fit, ifelse(df_bar$fit < 0, 0, df_bar$fit))
df_bar$num2 <-  ifelse(df_bar$num > opt$gnum, opt$gnum, ifelse(df_bar$num < 0, 0, df_bar$num))

# barplot
p2 <- NULL
p2 <- ggplot(df_bar, aes(x = as.factor(ALBST), y = num)) +
  geom_col(aes(fill = num2), width = 1) +
  theme_light()+
  scale_fill_gradient(low = "white", high = "#336699", limits = c(0, opt$gnum)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = 6,angle = 90, hjust = 1))+
  geom_line(aes(y = fit, colour = fit2, group = "ALBST"), size = 1.7)+
  scale_colour_gradient(low = "red", high = "#FFFF33", limits = c(0, opt$fit) )+
  coord_cartesian( ylim=c(opt$miny, NA))+
  ggtitle("expression_ratio")+
  theme(aspect.ratio = 0.15)+
  labs(x = NULL, y = "expression_ratio" )

# dotplot
meta <- obj@meta.data
p3 <- ggplot(data = meta, aes(x = ALBST, y = Width))+
  geom_point(colour = "black", size = 0.0005) +
  theme_light()+
  theme(aspect.ratio = 0.03, panel.grid = element_blank(), panel.border = element_blank(), axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())+
  labs(x = "", y = "")

pdf(file=paste(opt$outdir,"/",opt$prefix,"_spatial_heatmap.pdf", sep=""), width = w, height = h)
grid.arrange(p2, p1_grob, p3, ncol = 1)
dev.off()

png(filename=paste(opt$outdir,"/",opt$prefix,"_spatial_heatmap.png", sep=""), width = w*300, height = h*300, res=300, units="px")
grid.arrange(p2, p1_grob, p3, ncol = 1)
dev.off()








