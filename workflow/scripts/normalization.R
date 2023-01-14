# Normalization of gene expression data
# By João Sequeira
# Jan 2019

library("optparse")

option_list = list(
    make_option(c("-c", "--counts"), type="character", default=NULL,
                help="Table with abundance information.", metavar="character"),
    make_option(c("-m", "--method"), type="character", default="auto",
                help="Normalization method to apply (TMM, RLE or VSN). If 'auto', will determine from filename",
                metavar="character"),
      make_option(c("-o", "--output"), type="character", default=NULL,
              help="Output filename of normalized counts.", metavar="character")
);
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

print("Reading data to normalize.")
df <- read.table(opt$counts, header=TRUE, sep="\t", row.names=1)

# TMM or RLE normalization -> for RNA-Seq
if(opt$method == "TMM" || opt$method == "RLE") {
  print(paste("Performing", if (opt$method == "TMM") {"Trimmed Mean of M-values" } else {"Relative Log Expression"},
              "normalization.", sep=' '))
  library("edgeR")
  df[is.na(df)] <- 0
  # Remove
  factors = calcNormFactors(df, method = opt$method)
  write.table(factors, file = paste0(dirname(opt$counts), "/norm_factors.txt"), sep='\n',
              row.names = FALSE, col.names = FALSE)
  df[,1:ncol(df)] <- mapply("*", df[,1:ncol(df)], factors)

# Variance Stabilizing Normalization -> for proteomics
} else if(opt$method == "VSN") {
  df = as.matrix(df)              # ExpressionSet requires a matrix
  print("Performing Variance Stabilizing Normalization.")
  library("vsn")          # justvsn
  library("pcaMethods")   # llsImpute
  df[df==0] = NA                  # Replace 0s with NAs for next functions
  norm = exprs(justvsn(ExpressionSet(df)))

  # Local Least Squares Imputation
  print("Performing Local Least Squares Imputation.")
  allVariables = FALSE
  if (sum(complete.cases(norm)) / nrow(norm) < 0.5) {allVariables = TRUE}
  imputed = llsImpute(t(norm), correlation = "pearson", allVariables = allVariables)
  df = t(completeObs(imputed))

} else {
  stop("Error: normalization method must be either TMM, RLE or VSN")
}

print('Writting normalized results.')
write.table(df, file=opt$output, sep = "\t", row.names = TRUE,  col.names = TRUE)
