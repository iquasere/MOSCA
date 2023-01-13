# Differential expression analysis
# By João Sequeira
# Sep 2017

# problems with libreadline.so.6 might be solved with cd /lib/x86_64-linux-gnu/; sudo ln -s libreadline.so.7.0 libreadline.so.6

packages <- c("optparse", "DESeq2", "pheatmap", "RColorBrewer", "ROTS")

for (package in packages){
  eval(bquote(library(.(package))))
}

option_list <- list(
  make_option(c("-r", "--counts"), type="character", default=NULL, metavar="character",
              help="The expression matrix"),
  make_option(c("-c", "--conditions"), type="list", metavar="character", default=NULL,
              help="The conditions to define duplicates (e.g. 'c1,c1,c2,c2')"),
  make_option(c("-m", "--method"), default="differential", type="character", metavar="character",
              help="Method for ordering rows in protein expression heatmap
                [differential/abundance]"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="Output directory", metavar="character"),
  make_option(c("-f", "--foldchange"), type="numeric", default=2,
              help="Minimum fold change for two-tailed hypothesis of differential expression significance"),
  make_option(c("-d", "--data-type"), type="character", default="rna_seq", metavar="character",
              help="Data type [rna_seq/proteomics]"),
  make_option("--fdr", type="numeric", default=0.05, help="False discovery rate for differential analysis"))

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

paste("Counts:", opt$counts, sep=' ')
paste("Conditions:", opt$conditions, sep=' ')
paste("Method:", opt$method, sep=' ')
paste("Output:", opt$output, sep=' ')
paste("Minimum fold change:", opt$foldchange, sep=' ')

# Input parsing
opt$conditions <- strsplit(opt$conditions, ",")[[1]]
total <- read.table(opt$counts, h=T, row.names=1, sep = '\t')
conditions <- factor(opt$conditions)

# RNA-Seq differential analysis
if(opt$data_type == "rna_seq") {
  total <- total[ rowSums(total) > 1, ]
  cd <- data.frame(opt$conditions)
  colnames(cd)[1] <- "condition"
  rownames(cd) <- colnames(total)

  dds <- DESeqDataSetFromMatrix(countData = total, colData = cd, design = ~condition)
  dds <- DESeq(dds)
  res <- results(dds, lfcThreshold = log2(opt$foldchange))
  data <- counts(estimateSizeFactors(dds), normalized=TRUE)
  write.table(data, file=paste(file=opt$output, "normalized_counts.tsv", sep = '/'), sep='\t', col.names = NA, quote=FALSE)

  # Bland–Altman plot
  jpeg(paste(opt$output, "ma.jpeg", sep = '/'))
  plotMA(res, main="DESeq2", ylim=c(-2,2))
  dev.off()

  # Normalized counts
  jpeg(paste(opt$output, "counts.jpeg", sep = '/'))
  plotCounts(dds, gene=which.min(res$padj), intgroup="condition")
  dev.off()

  # Protein expressions differential analysis
  if(identical(opt$method, "differential")) {
    resOrdered <- res[order(res$padj),]
  } else {
    resOrdered <- res[order(-res$baseMean),]}
  write.table(as.data.frame(resOrdered), file=paste(file=opt$output, "condition_treated_results.tsv", sep = '/'),
              sep='\t', col.names = NA, quote=FALSE)
  vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
  select <- rownames(head(resOrdered,20))
  vsd.counts <- assay(vsd)[select,]
  jpeg(paste(opt$output, "gene_expression.jpeg", sep = '/'))
  pheatmap(vsd.counts)
  dev.off()

  # Sample expressions differential analysis
  sampleDists <- dist(t(assay(vsd)))
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- colnames(total)
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

# Proteomics differential analysis
} else if(opt$data_type == "proteomics") {
  # Reproducibility-Optimized Test Statistics
  de = ROTS(data=total, groups=opt$conditions, progress=TRUE)

  de_results <- cbind(de$logfc, de$pvalue)
  colnames(de_results) <- c("log2FoldChange", "pvalue")
  write.table(de_results, file=paste(opt$output, "condition_treated_results.tsv", sep = '/'), sep="\t", quote=FALSE)

  # Write parameters of DE model
  for (col in c("B", "a1", "a2", "k", "R", "Z")){
    write(sprintf("%s%s: %s", col, strrep(" ", 3-nchar(col)), de[col]), file= paste0(opt$output, '/report.txt'),
          append=TRUE)
  }

  #summary(de, fdr = 0.05)
  for (ptype in c('volcano', 'heatmap', 'ma', 'reproducibility', 'pvalue', 'pca')){
    jpeg(paste0(opt$output, "/", ptype, ".jpeg"))
    plot(de, type=ptype, fdr=opt$fdr)
    dev.off()
  }
}