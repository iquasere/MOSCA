# Title     : proteomics processing
# Objective : Normalize proteomics quantification and input missing values
# Created by: Joao Sequeira
# Created on: 08/12/2020

packages <- c("optparse",
              "vsn",
              "pcaMethods", # llsImpute
              "ROTS",
              "progress")   # progressbar

for (package in packages){
    eval(bquote(library(.(package))))
}

option_list = list(
  make_option(c("-m", "--matrix"), help="The quantification matrix (TSV format)", metavar="FILENAME"),
  make_option(c("-k", "--cluster-size"), type="list", metavar="N OF PROTEINS", default=10,
                help="Number of proteins to use for regression"),
  make_option(c("-o", "--output"), help="Output file", metavar="OUT FILENAME"),
  make_option(c("-c", "--conditions"), metavar="COMMA-SEPARATED LIST",
                help="The conditions to define replicates for differential analysis (e.g. 'c1,c1,c2,c2')"),
  make_option("--fdr", type="numeric", default=0.05, help="False discovery rate for differential analysis"))

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
print("Parameters inputted.")
opt$conditions <- strsplit(opt$conditions, ",")[[1]]

remove_nas_by_condition <- function (data, conditions) {
  vals = conditions[!duplicated(conditions)]
  cols = list()
  for (val in vals) {
    cols = append(cols, list(which(conditions == val)))
  }
  rows_to_keep = rownames(data)
  for (col_group in cols) {
    pb <- progress_bar$new(total = length(rows_to_keep), format = "Removing rows with < 2 non NA per sample: :current/:total [:bar] :percent eta: :eta")
    for (row in rows_to_keep) {
      if (sum(!is.na(data[row, colnames(data)[col_group]])) < 2) {
        rows_to_keep = rows_to_keep[rows_to_keep != row]
      }
      pb$tick()
    }}
  return(data[rows_to_keep,])
}

pdata = as.matrix(read.table(opt$matrix, sep='\t', h=T, row.names=1))
pdata[pdata==0] = NA    # maybe add 1 in the future, if LLS imputation keeps failling?
print("Matrix read.")

# variance stabilization normalization
norm = justvsn(ExpressionSet(pdata))
print("Normalization done.")

filtered = remove_nas_by_condition(exprs(norm), opt$conditions)  # remove rows with more than 1 NA per sample
#normalized = normalized[rowSums(is.na(normalized)) < 2, , drop=FALSE]   # love R's insanity
print("Filtered rows with too many NAs.")

# local least squares imputation - TODO - check when this works
#allVariables = TRUE
# using 'nniRes' function is advised in documentation, but it does not allow to set 'allVariables' to TRUE
# not using 'allVariables' will result in too many missing values breaking the script
#imputed = llsImpute(t(exprs(norm)), correlation = "pearson", allVariables = allVariables)
#observations = t(completeObs(imputed))
#print("Imputation done.")

# export results
write.table(filtered, file = paste0(opt$output, '/quantification.tsv'), sep='\t', row.names = TRUE, col.names = TRUE)
print("Results exported.")

# Reproducibility-Optimized Test Statistic
de = ROTS(data=filtered, groups=opt$conditions, progress=TRUE)
for (col in c("B", "a1", "a2", "k", "R", "Z")){
  write(sprintf("%s%s: %s", col, strrep(" ", 3-nchar(col)), de[col]), file= paste0(opt$output, '/report.txt'),
        append=TRUE)
}
print("Differential analysis done.")

#summary(de, fdr = 0.05)
for (ptype in c('volcano', 'heatmap', 'ma', 'reproducibility', 'pvalue', 'pca')){
  jpeg(paste0(opt$output, "/", ptype, ".jpeg"))
  plot(de, type=ptype, fdr=opt$fdr)
  dev.off()
}
