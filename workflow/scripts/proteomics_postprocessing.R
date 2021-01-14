# Title     : proteomics processing
# Objective : Normalize proteomics quantification and input missing values
# Created by: Joao Sequeira
# Created on: 08/12/2020

packages <- c("optparse", "vsn", "pcaMethods", "ROTS")

for (package in packages){
    eval(bquote(library(.(package))))
}

option_list = list(
    make_option(c("-qm", "--quantification-matrix"), help="The quantification matrix (TSV format)", metavar="character"),
    make_option(c("-k", "--cluster-size"), type="list", metavar="character", default=10,
                help="Number of proteins to use for regression"),
    make_option(c("-o", "--output"), type="character", help="Output file", metavar="character"),
    make_option(c("-c", "--conditions"), type="list",
                help="The conditions to define replicates for differential analysis (e.g. 'c1,c1,c2,c2')"));

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

opt$conditions <- strsplit(opt$conditions, ",")[[1]]

pdata = as.matrix(read.table(opt$quantification_matrix, sep='\t', h=T, row.names=1))
pdata[pdata==0] = NA

# variance stabilization normalization
norm = justvsn(ExpressionSet(pdata))

# local least squares imputation
imputed = llsImpute(t(exprs(norm)))
observations = t(completeObs(imputed))

# export results
write.table(new_data, file = opt$output, sep='\t', row.names = TRUE, col.names = TRUE)

# Reproducibility-Optimized Test Statistic
de = ROTS(data=observations, groups=conditions, progress=TRUE)

for (col in c("B", "a1", "a2", "k", "R", "Z")){
  write(sprintf("%s%s: %s", col, strrep(" ", 3-nchar(col)), de[col]), file='report.txt', append=TRUE)
}

