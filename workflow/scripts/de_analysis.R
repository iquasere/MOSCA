# Differential expression analysis
# By João Sequeira
# Sep 2017

# problems with libreadline.so.6 might be solved with cd /lib/x86_64-linux-gnu/; sudo ln -s libreadline.so.7.0 libreadline.so.6
method <- 'differential' # snakemake@params$method TODO - think about reintegrating the method by baseMean

paste0("Counts: ", snakemake@input[[1]])
paste0("Conditions: ", snakemake@params$conditions)
paste0("Method: ", method)
paste0("Output: ", snakemake@params$output)
paste0("Minimum fold change: ", snakemake@params$foldchange)

# Input parsing
conditions <- strsplit(snakemake@params$conditions, ",")[[1]]
total <- read.table(snakemake@input[[1]], h=T, row.names=1, sep = '\t')

# RNA-Seq differential analysis
if(snakemake@params$datatype == "rna_seq") {
  for (package in c("DESeq2", "pheatmap", "RColorBrewer")){
    eval(bquote(library(.(package))))
  }
  conditions <- factor(conditions)
  total[is.na(total)] <- 0
  total <- total[ rowSums(total) > 1, ]
  cd <- data.frame(conditions)
  colnames(cd)[1] <- "condition"
  rownames(cd) <- colnames(total)

  dds <- DESeqDataSetFromMatrix(countData = total, colData = cd, design = ~condition)
  dds <- DESeq(dds)
  res <- results(dds, lfcThreshold = log2(snakemake@params$foldchange))
  data <- counts(estimateSizeFactors(dds), normalized=TRUE)
  write.table(
    data, file=paste0(file=snakemake@params$output, "/normalized_counts.tsv"), sep='\t', col.names = NA, quote=FALSE)

  # Bland–Altman plot
  jpeg(paste0(snakemake@params$output, "/ma.jpeg"))
  plotMA(res, main="DESeq2", ylim=c(-2,2))
  dev.off()

  # Normalized counts
  jpeg(paste0(snakemake@params$output, "/counts.jpeg"))
  plotCounts(dds, gene=which.min(res$padj), intgroup="condition")
  dev.off()

  # Protein expressions differential analysis
  if(identical(method, "differential")) {
    resOrdered <- res[order(res$padj),]
  } else {
    resOrdered <- res[order(-res$baseMean),]}
  write.table(as.data.frame(resOrdered), file=paste0(
    snakemake@params$output, "/condition_treated_results.tsv"), sep='\t', col.names = NA, quote=FALSE)
  vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
  select <- rownames(head(resOrdered,30))
  vsd.counts <- assay(vsd)[select,]
  jpeg(paste0(snakemake@params$output, "/gene_expression.jpeg"))
  pheatmap(vsd.counts)
  dev.off()

  # Sample expressions differential analysis
  sampleDists <- dist(t(assay(vsd)))
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- colnames(total)
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  jpeg(paste0(snakemake@params$output, "/sample_distances.jpeg"))
  pheatmap(
    sampleDistMatrix, clustering_distance_rows = sampleDists, clustering_distance_cols = sampleDists, col = colors)
  dev.off()

  # Principal components analysis
  jpeg(paste0(snakemake@params$output, "/pca.jpeg"))
  plotPCA(vsd, intgroup = c("condition"))
  dev.off()

# Proteomics differential analysis
} else if(snakemake@params$datatype == "proteomics") {
  library("ROTS")
  library("tibble")
  # Reproducibility-Optimized Test Statistics
  de <- ROTS(data=total, groups=conditions, progress=TRUE)
  de_res <- cbind(de$logfc, de$pvalue, de$FDR)
  colnames(de_res) <- c("log2FoldChange", "pvalue", "FDR")
  sorted_res <- de_res[order(de_res[, "pvalue"]), ]
  write.table(
    tibble::rownames_to_column(as.data.frame(sorted_res), "IDs"),     # add rownames as column
    file=paste0(snakemake@params$output, "/condition_treated_results.tsv"),
    sep="\t", quote=FALSE, row.names=FALSE)

  # Write parameters of DE model
  for (col in c("B", "a1", "a2", "k", "R", "Z")){
    write(sprintf("%s%s: %s", col, strrep(" ", 3-nchar(col)), de[col]),
          file=paste0(snakemake@params$output, '/report.txt'), append=TRUE)
  }

  sorted_data <- head(as.matrix(total[order(de$pvalue), ]), 30)
  #summary(de, fdr = 0.05)
  for (ptype in c('volcano', 'heatmap', 'ma', 'reproducibility', 'pvalue', 'pca')) {
    print(paste0("Plotting ", ptype, " plot."))
    jpeg_file <- paste0(ptype, ".jpeg")

    tryCatch({
      if (ptype == 'heatmap') {
        pheatmap(sorted_data, Rowv = NA, Colv = NA, col = colorRampPalette(c("blue", "white", "red"))(256), legend = TRUE)
      } else {
        jpeg(jpeg_file, res = 300)
        plot(de, type = ptype, fdr = 0.05)
        dev.off()
      }
    }, error = function(e) {
      tryCatch({
        if (ptype == 'heatmap') {
          pheatmap(sorted_data, Rowv = NA, Colv = NA, col = colorRampPalette(c("blue", "white", "red"))(256), legend = TRUE)
        } else {
          jpeg(jpeg_file)     # If plot generation fails due to margin size, try without specifying res
          plot(de, type = ptype, fdr = 0.05)
          dev.off()
        }
      }, error = function(inner_e) {
        cat("Error in plotting", ptype, "plot:", conditionMessage(inner_e), "\n")
      })
    })
  }
}

print("DE analysis finished.")