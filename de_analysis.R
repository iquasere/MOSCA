# Differential expression analysis
# By João Sequeira
# Sep 2017

library("optparse")
library("DESeq2")
library("pheatmap")
library("RColorBrewer")

option_list = list(
    make_option(c("-r", "--readcounts"), type="character", default=NULL, 
                help="The expression matrix", metavar="character"),
    make_option(c("-c", "--conditions"), type="list", metavar="character",
                help="The conditions to define duplicates (e.g. 'c1 c1 c2 c2')",
                default=NULL),
    make_option(c("-m,", "--method"), default="differential", type="character", 
                help="Method for ordering rows in protein expression heatmap 
                [differential/abundance]", metavar="character",),
    make_option(c("-o", "--output"), type="character", default=NULL, 
                help="Output directory", metavar="character"));
    
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

total <- read.table(opt$readcounts, h=T, row.names=1)
condition <- factor(opt$conditions)
total <- total[ rowSums(total) > 1, ]
cd = data.frame(opt$conditions)
colnames(cd)[1]="condition"
rownames(cd)=colnames(total)
dds <- DESeqDataSetFromMatrix(countData = total, colData = cd, design = ~condition)
dds <- DESeq(dds)
res <- results(dds)
resOrdered <- ifelse(opt$method == "differential", res[order(res$padj),], 
                                                res[order(-res$baseMean),])

jpeg(opt$output + "/ma.jpeg")
plotMA(res, main="DESeq2", ylim=c(-2,2))
dev.off()

jpeg(opt$output + "/counts.jpeg")
plotCounts(dds, gene=which.min(res$padj), intgroup="condition")
dev.off()

write.csv(as.data.frame(resOrdered), file=opt$output + "/condition_treated_results.csv")

vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
select=rownames(head(resOrdered,20))
vsd.counts = assay(vsd)[select,]
jpeg(opt$output + "/gene_expression.jpeg")
pheatmap(vsd.counts)
dev.off()

sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- dds$condition
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
jpeg(opt$output + "/sample_distances.jpeg")
pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors)
dev.off()

jpeg(opt$output + "/pca.jpeg")
plotPCA(vsd, intgroup=c("condition"))
dev.off()