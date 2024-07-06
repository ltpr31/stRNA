library(Seurat)
library(monocle)
library(ggpubr)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library("argparse")
parser <- ArgumentParser(description='cooperation with monocle')
parser$add_argument("-i", "--rds", type="character",required=T)
parser$add_argument( "-n", "--heatmap.gene", type="integer", default=20)
parser$add_argument( "-c", "--heatmap.clusters", type="integer", default=6)
parser$add_argument( "-H", "--height", type="double", default=7)
parser$add_argument("-W", "--width", type="double", default=6)
parser$add_argument( "-o", "--outdir", type="character", default=getwd())
parser$add_argument("-p", "--prefix", type="character", default="demo")

opt <- parser$parse_args()
if( !file.exists(opt$outdir) ){
	if( !dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE) ){
		stop(paste("failed to create path: ",opt$outdir,sep=""))
	}
}


qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))) 

options(future.globals.maxSize = 500 * 1024^3)
stRNA <- readRDS(opt$rds)
DefaultAssay(stRNA) <- "Spatial"
message(DefaultAssay(stRNA))
assay="Spatial"
head(stRNA@meta.data)
expr_matrix <- stRNA@assays[[assay]]@counts
p_data <- stRNA@meta.data
head(p_data)
f_data <- data.frame(gene_short_name = row.names(stRNA),row.names = row.names(stRNA))
pd <- new('AnnotatedDataFrame', data = p_data)
fd <- new('AnnotatedDataFrame', data = f_data)

cds <- newCellDataSet(expr_matrix,
		phenoData = pd,
		featureData = fd,
		lowerDetectionLimit = 0.5,
		expressionFamily = negbinomial.size())

cds<-estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

######### calculation of pseudotime

disp_table <- dispersionTable(cds)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
cds <- setOrderingFilter(cds, disp.genes)
high_var_genes_info <- disp_table[disp_table$gene_id %in% disp.genes, ]
write.table(high_var_genes_info, file=paste(opt$outdir,"/",opt$prefix,".monocle_highvarible.Gene.info.tsv",sep=""), row.names = F, col.names = T, quote = F)
p1<-plot_ordering_genes(cds)
pdf(file=paste(opt$outdir,"/",opt$prefix,".selected_order_gene.pdf",sep=""), height=opt$height, width=opt$width)
print(p1)
dev.off()
png(filename=paste(opt$outdir,"/",opt$prefix,".selected_order_gene.png",sep=""), height=opt$height*300, width=opt$width*300, res=300, units="px")
print(p1)
dev.off()
rm(stRNA)
cds <- reduceDimension(cds, max_components = 2,
		method = 'DDRTree')
cds <- orderCells(cds)

###################dimensionality reduction with ALBST and Width：
reduction <- p_data[,c("ALBST","Width")]
reduction <- t(reduction)
cds@reducedDimS <- reduction
cds@reducedDimK <- reduction
saveRDS(cds, file = paste(opt$outdir,"/",opt$prefix,".spatial_trajectory.rds",sep=""))
df<-data.frame(barcode=row.names(pData(cds)),pData(cds))
write.table(df,file=paste(opt$outdir,"/",opt$prefix,".Pseudotime.metadata.tsv",sep=""),sep="\t",quote=F,row.names=F,col.names = T)
p=plot_cell_trajectory(cds, color_by = "Pseudotime", cell_size = 0.3)
pdf(file=paste(opt$outdir,"/",opt$prefix,".Pseudotime.pdf",sep=""), height=opt$height*2/7, width=opt$width)
p
dev.off()
png(filename=paste(opt$outdir,"/",opt$prefix,".Pseudotime.png",sep=""), height=opt$height*300*2/7, width=opt$width*300, res=300, units="px")
p
dev.off()
p=plot_cell_trajectory(cds, color_by = "State", cell_size = 0.3)+scale_colour_manual(values =  col_vector)
pdf(file="State.pdf", height=opt$height*2/7, width=opt$width)
p
dev.off()
png(filename=paste(opt$outdir,"/",opt$prefix,".State.png",sep=""), height=opt$height*300*2/7, width=opt$width*300, res=300, units="px")
p
dev.off()


########## replace the pseudotime with ALBST：
pData(cds)$Pseudotime <- pData(cds)$ALBST


########## identify ALBST related genes:
Time_diff <- differentialGeneTest(cds[disp.genes,], cores = 1,fullModelFormulaStr = "~sm.ns(Pseudotime)")
Time_diff=arrange(Time_diff,qval)
write.table(data.frame(ID=rownames(Time_diff),Time_diff),file=paste(opt$outdir,"/",opt$prefix,".albst_diff_all.tsv",sep=""),sep="\t",quote=F,row.names=F,col.names = T)
Time_genes <- top_n(Time_diff, n = opt$heatmap.gene, desc(qval)) %>% pull(gene_short_name) %>% as.character()
pdf(file=paste(opt$outdir,"/",opt$prefix,".topALBST_diff_genes_heatmap.pdf",sep=""), height=opt$height*2, width=opt$width)
p=plot_pseudotime_heatmap(cds[Time_genes,],
		num_clusters = opt$heatmap.clusters,
		cores = 1,
		show_rownames = T, return_heatmap=T)
dev.off()
png(filename=paste(opt$outdir,"/",opt$prefix,".topALBST_diff_genes_heatmap.png",sep=""), height=opt$height*300*2, width=opt$width*300, res=300, units="px")
plot_pseudotime_heatmap(cds[Time_genes,],
		num_clusters = opt$heatmap.clusters,
		cores = 1,
		show_rownames = T)
dev.off()
clusters <- cutree(p$tree_row, k = opt$heatmap.clusters)
clustering <- data.frame(clusters)
clustering[,1] <- as.character(clustering[,1])
colnames(clustering) <- "Gene_Clusters"
write.table(data.frame(ID=rownames(clustering),clustering),file=paste(opt$outdir,"/",opt$prefix,".topALBST_diff_genes_heatmap_cluster.tsv",sep=""),sep="\t",quote=F,row.names=F,col.names = T)




