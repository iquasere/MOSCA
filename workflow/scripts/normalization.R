# Normalization of gene expression data
# By João Sequeira
# Jan 2019

normalize_gene_expression <- function(df, output_file, method) {
  # TMM or RLE normalization -> for RNA-Seq
  if (method == "TMM" || method == "RLE") {
    print(paste("Performing", if (method == "TMM") {"Trimmed Mean of M-values"} else {"Relative Log Expression"}, "normalization.", sep=' '))
    library("edgeR")
    df[is.na(df)] <- 0
    factors <- calcNormFactors(df, method=method)
    write.table(factors, file=paste0(dirname(output_file), "/norm_factors.txt"), sep='\n', row.names=FALSE, col.names=FALSE)
    df[, 1:ncol(df)] <- mapply("*", df[, 1:ncol(df)], factors)
  } else if (method == "VSN") {
    df <- as.matrix(df)
    print("Performing Variance Stabilizing Normalization.")
    library("vsn")
    library("pcaMethods")
    df[df == 0] <- NA
    norm <- exprs(justvsn(ExpressionSet(df)))
    print("Performing Local Least Squares Imputation.")
    allVariables <- FALSE
    if (sum(complete.cases(norm)) / nrow(norm) < 0.5) {
      allVariables <- TRUE
    }
    imputed <- llsImpute(t(norm), correlation = "pearson", allVariables = allVariables)
    df <- t(completeObs(imputed))
  } else {
    stop("Error: normalization method must be either TMM, RLE, or VSN")
  }
  print('Writing normalized results.')
  write.table(df, file = output_file, sep = "\t", row.names = TRUE, col.names = TRUE)
}

print("Reading data to normalize.")
df1 <- read.table(snakemake@input[[1]], header=TRUE, sep="\t", row.names=1)
normalize_gene_expression(df1, snakemake@output[[1]], snakemake@params$method)
if (length(snakemake@input) > 1) {
  print("Reading more data to normalize.")
  df2 <- read.table(snakemake@input[[2]], header=TRUE, sep="\t", row.names=1)
  # check if mt is in the name (snakemake@input[[2]]) and mp is not
  if (grepl("mt", snakemake@input[[2]]) && !grepl("mp", snakemake@input[[2]])) {
      normalize_gene_expression(df2, snakemake@output[[2]], snakemake@params$method)
  } else {
      normalize_gene_expression(df2, snakemake@output[[2]], "VSN")
  }
}
