# Title     : proteomics processing
# Objective : Normalize proteomics quantification and input missing values
# Created by: Joao Sequeira
# Created on: 08/12/2020

packages <- c("optparse", "vsn", "pcaMethods", "ROTS")

for (package in packages){
    eval(bquote(library(.(package))))
}

option_list = list(
    make_option(c("-m", "--matrix"), help="The quantification matrix (TSV format)", metavar="FILENAME"),
    make_option(c("-k", "--cluster-size"), type="list", metavar="N OF PROTEINS", default=10,
                help="Number of proteins to use for regression"),
    make_option(c("-o", "--output"), help="Output file", metavar="OUT FILENAME"),
    make_option(c("-c", "--conditions"), metavar="COMMA-SEPARATED LIST",
                help="The conditions to define replicates for differential analysis (e.g. 'c1,c1,c2,c2')"));

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
print("Parameters inputted.")

opt$conditions <- strsplit(opt$conditions, ",")[[1]]

pdata = as.matrix(read.table(opt$matrix, sep='\t', h=T, row.names=1))
pdata[pdata==0] = NA
print("Matrix read.")

# variance stabilization normalization
norm = justvsn(ExpressionSet(pdata))
print("Normalization done.")

# local least squares imputation
imputed = llsImpute(t(exprs(norm)))
observations = t(completeObs(imputed))
print("Imputation done.")

# export results
write.table(new_data, file = opt$output, sep='\t', row.names = TRUE, col.names = TRUE)
print("Results exported.")

# Reproducibility-Optimized Test Statistic
de = ROTS(data=observations, groups=conditions, progress=TRUE)
print("Differential analysis done.")

for (col in c("B", "a1", "a2", "k", "R", "Z")){
  write(sprintf("%s%s: %s", col, strrep(" ", 3-nchar(col)), de[col]), file='report.txt', append=TRUE)
}
