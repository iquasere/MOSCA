# Differential expression analysis
# By João Sequeira
# Sep 2017

library("optparse")

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

paste("Readcounts:", opt$readcounts, sep=' ')
paste("Conditions:", opt$conditions, sep=' ')
paste("Method:", opt$method, sep=' ')
paste("Output:", opt$output, sep=' ')

opt$conditions <- strsplit(opt$conditions, "[[:space:]]")[[1]]

total <- read.table(opt$readcounts, h=T, row.names=1, sep = '\t')
condition <- factor(opt$conditions)
total <- total[ rowSums(total) > 1, ]
cd = data.frame(opt$conditions)
colnames(cd)[1]="condition"
rownames(cd)=colnames(total)
library("DESeq2")

dds <- DESeqDataSetFromMatrix(countData = total, colData = cd, design = ~condition)
dds <- DESeq(dds)
res <- results(dds)

# Bland–Altman plot
jpeg(paste(opt$output, "ma.jpeg", sep = '/'))
plotMA(res, main="DESeq2", ylim=c(-2,2))
dev.off()

# Normalized counts
jpeg(paste(opt$output, "counts.jpeg", sep = '/'))
plotCounts(dds, gene=which.min(res$padj), intgroup="condition")
dev.off()

library("pheatmap")
library("RColorBrewer")

# Protein expressions differential analysis
if(identical(opt$method, "differential")) {
    resOrdered <- res[order(res$padj),]
} else {
    resOrdered <- res[order(-res$baseMean),]}
write.csv(as.data.frame(resOrdered), paste(file=opt$output, "condition_treated_results.csv", sep = '/'))
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
select=rownames(head(resOrdered,20))
vsd.counts = assay(vsd)[select,]
jpeg(paste(opt$output, "gene_expression.jpeg", sep = '/'))
pheatmap(vsd.counts)
dev.off()

# Sample expressions differential analysis
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- dds$condition
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
jpeg(paste(opt$output, "sample_distances.jpeg", sep = '/'))
pheatmap(sampleDistMatrix, clustering_distance_rows = sampleDists, 
         clustering_distance_cols = sampleDists, col = colors)
dev.off()

# Principal components analysis
jpeg(paste(opt$output, "pca.jpeg", sep = '/'))
plotPCA(vsd, intgroup = c("condition"))
dev.off()