# TMM normalization
# By João Sequeira
# Jan 2019


library("optparse")
library("edgeR")

option_list = list(
    make_option(c("-t", "--table"), type="character", default=NULL, 
                help="table with abundance information", metavar="character"),
    make_option(c("-o", "--output"), type="character", default="factors.txt", 
                help="output file name [default= %default]", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

df <- read.table(opt$table, header = TRUE, sep = "\t")

if(sapply(df, class)[1] != "numeric" & sapply(df, class)[1] != "integer"){
        df[colnames(df)[1]] <- NULL
        }

normalization_factors = calcNormFactors(df)

write.table(normalization_factors, file = opt$output, sep='\t', row.names = FALSE, col.names = FALSE)
