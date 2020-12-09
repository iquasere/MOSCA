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
    make_option(c("-c", "--conditions"), type="list", metavar="character",
                help="The conditions to define duplicates (e.g. 'c1,c1,c2,c2')",
                default=NULL),
    make_option(c("-m", "--method"), default="differential", type="character",
                help="Method for ordering rows in protein expression heatmap
                [differential/abundance]", metavar="character",),
    make_option(c("-o", "--output"), type="character", default=NULL,
                help="Output directory", metavar="character"));

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

pdata = as.matrix(read.table(opt$quantification_matrix, sep='\t', h=T, row.names=1))
pdata[pdata==0] = NA
pset = ExpressionSet(pdata)

norm = justvsn(pset)
imputed = llsImpute(pset)
