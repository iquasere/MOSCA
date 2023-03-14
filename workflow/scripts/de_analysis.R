# Differential expression analysis
# By João Sequeira
# Sep 2017

# problems with libreadline.so.6 might be solved with cd /lib/x86_64-linux-gnu/; sudo ln -s libreadline.so.7.0 libreadline.so.6

paste("Counts:", snakemake@input[[1]], sep=' ')
paste("Conditions:", snakemake@params$conditions, sep=' ')
paste("Method:", snakemake@params$method, sep=' ')
paste("Output:", snakemake@params$output, sep=' ')
paste("Minimum fold change:", snakemake@params$foldchange, sep=' ')

# Input parsing
conditions <- strsplit(snakemake@params$conditions, ",")[[1]]
total <- read.table(snakemake@input[[1]], h=T, row.names=1, sep = '\t')
conditions <- factor(snakemake@params$conditions)

# RNA-Seq differential analysis
if(snakemake@params$datatype == "rna_seq") {
  for (package in c("DESeq2", "pheatmap", "RColorBrewer")){
    eval(bquote(library(.(package))))
  }
  total[is.na(total)] <- 0
  total <- total[ rowSums(total) > 1, ]
  cd <- data.frame(conditions)
  colnames(cd)[1] <- "condition"
  rownames(cd) <- colnames(total)

  dds <- DESeqDataSetFromMatrix(countData = total, colData = cd, design = ~condition)
  dds <- DESeq(dds)
  res <- results(dds, lfcThreshold = log2(snakemake@params$foldchange))
  data <- counts(estimateSizeFactors(dds), normalized=TRUE)
  write.table(data, file=paste(file=snakemake@params$output, "normalized_counts.tsv", sep = '/'),
              sep='\t', col.names = NA, quote=FALSE)

  # Bland–Altman plot
  jpeg(paste(snakemake@params$output, "ma.jpeg", sep = '/'))
  plotMA(res, main="DESeq2", ylim=c(-2,2))
  dev.off()

  # Normalized counts
  jpeg(paste(snakemake@params$output, "counts.jpeg", sep = '/'))
  plotCounts(dds, gene=which.min(res$padj), intgroup="condition")
  dev.off()

  # Protein expressions differential analysis
  if(identical(snakemake@params$method, "differential")) {
    resOrdered <- res[order(res$padj),]
  } else {
    resOrdered <- res[order(-res$baseMean),]}
  write.table(as.data.frame(resOrdered), file=paste(
    file=snakemake@params$output, "condition_treated_results.tsv", sep = '/'), sep='\t', col.names = NA, quote=FALSE)
  vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
  select <- rownames(head(resOrdered,20))
  vsd.counts <- assay(vsd)[select,]
  jpeg(paste(snakemake@params$output, "gene_expression.jpeg", sep = '/'))
  pheatmap(vsd.counts)
  dev.off()

  # Sample expressions differential analysis
  sampleDists <- dist(t(assay(vsd)))
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- colnames(total)
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  jpeg(paste(snakemake@params$output, "sample_distances.jpeg", sep = '/'))
  pheatmap(sampleDistMatrix, clustering_distance_rows = sampleDists,
           clustering_distance_cols = sampleDists, col = colors)
  dev.off()

  # Principal components analysis
  jpeg(paste(snakemake@params$output, "pca.jpeg", sep = '/'))
  plotPCA(vsd, intgroup = c("condition"))
  dev.off()

# Proteomics differential analysis
} else if(snakemake@params$datatype == "proteomics") {
  library("ROTS")
  # Reproducibility-Optimized Test Statistics
  de = ROTS(data=total, groups=snakemake@params$conditions, progress=TRUE)

  de_results <- cbind(de$logfc, de$pvalue)
  colnames(de_results) <- c("log2FoldChange", "pvalue")
  write.table(
    de_results, file=paste(snakemake@params$output, "condition_treated_results.tsv", sep = '/'), sep="\t", quote=FALSE)

  # Write parameters of DE model
  for (col in c("B", "a1", "a2", "k", "R", "Z")){
    write(sprintf("%s%s: %s", col, strrep(" ", 3-nchar(col)), de[col]),
          file=paste0(snakemake@params$output, '/report.txt'), append=TRUE)
  }

  #summary(de, fdr = 0.05)
  for (ptype in c('volcano', 'heatmap', 'ma', 'reproducibility', 'pvalue', 'pca')){
    jpeg(paste0(snakemake@params$output, "/", ptype, ".jpeg"))
    plot(de, type=ptype, fdr=snakemake@params$fdr)
    dev.off()
  }
}