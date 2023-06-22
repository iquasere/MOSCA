# Normalization of gene expression data
# By João Sequeira
# Jan 2019

print("Reading data to normalize.")
df <- read.table(snakemake@input[[1]], header=TRUE, sep="\t", row.names=1)

# TMM or RLE normalization -> for RNA-Seq
if (snakemake@params$method == "TMM" || snakemake@params$method == "RLE") {
  print(paste(
    "Performing", if (snakemake@params$method == "TMM") {"Trimmed Mean of M-values" } else {"Relative Log Expression"},
    "normalization.", sep=' '))
  library("edgeR")
  df[is.na(df)] <- 0
  # Remove
  factors <- calcNormFactors(df, method=snakemake@params$method)
  write.table(factors, file=paste0(dirname(snakemake@output[[1]]), "/norm_factors.txt"), sep='\n',
              row.names=FALSE, col.names=FALSE)
  df[,1:ncol(df)] <- mapply("*", df[,1:ncol(df)], factors)

# Variance Stabilizing Normalization -> for proteomics
} else if(snakemake@params$method == "VSN") {
  df <- as.matrix(df)              # ExpressionSet requires a matrix
  print("Performing Variance Stabilizing Normalization.")
  library("vsn")          # justvsn
  library("pcaMethods")   # llsImpute
  df[df==0] <- NA                  # Replace 0s with NAs for next functions
  norm <- exprs(justvsn(ExpressionSet(df)))

  # Local Least Squares Imputation
  print("Performing Local Least Squares Imputation.")
  allVariables <- FALSE
  if (sum(complete.cases(norm)) / nrow(norm) < 0.5) {allVariables <- TRUE}
  imputed <- llsImpute(t(norm), correlation="pearson", allVariables=allVariables)
  df <- t(completeObs(imputed))

} else {
  stop("Error: normalization method must be either TMM, RLE or VSN")
}

print('Writting normalized results.')

write.table(df, file=snakemake@output[[1]], sep="\t", row.names=TRUE,  col.names=TRUE)
