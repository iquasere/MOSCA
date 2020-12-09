# Title     : proteomics processing
# Objective : Normalize proteomics quantification and input missing values
# Created by: Joao Sequeira
# Created on: 08/12/2020

packages <- c("optparse", "vsn", "pcaMethods")

for (package in packages){
    eval(bquote(library(.(package))))
    }

option_list = list(
    make_option(c("-qm", "--quantification-matrix"), type="character", default=NULL,
                help="The quantification matrix (TSV format)", metavar="character"),
    make_option(c("-k", "--cluster-size"), type="list", metavar="character",
                help="Number of proteins to use for regression", default=10),
    make_option(c("-o", "--output"), type="character", default=NULL,
                help="Output file", metavar="character"));

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

pdata = as.matrix(read.table(opt$quantification_matrix, sep='\t', h=T, row.names=1))
pdata[pdata==0] = NA

# variance stabilization normalization
norm = justvsn(ExpressionSet(pdata))

# local least squares imputation
imputed = llsImpute(t(exprs(norm)))
new_data = t(completeObs(imputed))

# export results
write.table(new_data, file = opt$output, sep='\t', row.names = TRUE, col.names = TRUE)