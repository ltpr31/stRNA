library(Seurat)
library(monocle)
library(ggpubr)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library("argparse")

parser <- ArgumentParser(description='monocle heatmap')

parser$add_argument("-i", "--cds", type="character",required=T)
parser$add_argument(  "--list", type="character",default=NULL)
parser$add_argument( "-n", "--heatmap.gene", type="integer", default=20000)
parser$add_argument( "-c", "--heatmap.clusters", type="integer", default=6)
parser$add_argument( "-H", "--height", type="double", default=7)
parser$add_argument("-W", "--width", type="double", default=6)
parser$add_argument( "-o", "--outdir", type="character", default=getwd())
parser$add_argument("-p", "--prefix", type="character", default="demo")
parser$add_argument("--no.cluster.row", action="store_true", default=FALSE)
opt <- parser$parse_args()

if( !file.exists(opt$outdir) ){
	if( !dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE) ){
		stop(paste("failed to create the path: ",opt$outdir,sep=""))
	}
}
options(future.globals.maxSize = 500 * 1024^3) 


cds <- readRDS(opt$cds)
head(pData(cds))

if (is.null(opt$list)){
disp_table <- dispersionTable(cds)
	disp.genes <- subset(disp_table, mean_expression >= 0.1)$gene_id
	high_var_genes_info <- disp_table[disp_table$gene_id %in% disp.genes, ]
	write.table(high_var_genes_info, file=paste(opt$outdir,"/",opt$prefix,".monocle_highvarible.Gene.info.tsv",sep=""), row.names = F, col.names = T, quote = F)

Time_diff <- differentialGeneTest(cds[disp.genes,], cores = 1,fullModelFormulaStr = "~sm.ns(Pseudotime)")
Time_diff=arrange(Time_diff,qval)
write.table(data.frame(ID=rownames(Time_diff),Time_diff),file=paste(opt$outdir,"/",opt$prefix,".Time_diff_all.tsv",sep=""),sep="\t",quote=F,row.names=F,col.names = T)

Time_genes <- top_n(Time_diff, n = opt$heatmap.gene, desc(qval)) %>% pull(gene_short_name) %>% as.character()
} else {
a <- read.table(opt$list, header = T)
a$ID <- gsub("_","-",a$ID)
a <- a$ID
Time_genes <- unique(a)
Time_genes <- Time_genes[Time_genes %in% rownames(cds)]
}

########## replace the pseudotime with ALBST：
pData(cds)$Pseudotime <- pData(cds)$ALBST
Time_genes <- Time_genes[Time_genes %in% rownames(cds)]
print(length(Time_genes))

if (opt$no.cluster.row){

pdf(file=paste(opt$outdir,"/",opt$prefix,".spatial.trajectory.heatmap.pdf",sep=""), height=opt$height*2, width=opt$width)
p=plot_pseudotime_heatmap(cds[Time_genes,],
		cores = 1,cluster_rows = FALSE,
		show_rownames = T, return_heatmap=T)
dev.off()

}else if (dim(cds[Time_genes,])[1] >= 2){
pdf(file=paste(opt$outdir,"/",opt$prefix,".spatial.trajectory.heatmap.pdf",sep=""), height=opt$height*2, width=opt$width)
p=plot_pseudotime_heatmap(cds[Time_genes,],
		num_clusters = opt$heatmap.clusters,
		cores = 1,cluster_rows = TRUE,
		show_rownames = T, return_heatmap=T)
dev.off()
}else{
pdf(file=paste(opt$outdir,"/",opt$prefix,".spatial.trajectory.heatmap.pdf",sep=""), height=opt$height*2, width=opt$width)
p=plot_pseudotime_heatmap(cds[Time_genes,],
		num_clusters = opt$heatmap.clusters,
		cores = 1,cluster_rows = FALSE,
		show_rownames = T, return_heatmap=T)
dev.off()
}


png(filename=paste(opt$outdir,"/",opt$prefix,".spatial.trajectory.heatmap.png",sep=""), height=opt$height*300*2, width=opt$width*300, res=300, units="px")
print(p)
dev.off()


clusters <- cutree(p$tree_row, k = opt$heatmap.clusters)
clustering <- data.frame(clusters)
clustering[,1] <- as.character(clustering[,1])
colnames(clustering) <- "Gene_Clusters"
write.table(data.frame(ID=rownames(clustering),clustering),file=paste(opt$outdir,"/",opt$prefix,".spatial.trajectory.cluster.tsv",sep=""),sep="\t",quote=F,row.names=F,col.names = T)






