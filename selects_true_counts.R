"Usage: marks_with_msp_bigger_than_0.R (--ms_insensitive_samples <mspI>) (--out1 <O1>) <input1>
-h --help    show this
--ms_insensitive_samples    <mspI>  Names of samples from Methylation Insensitive library.
--out1  name1   bed file with all the corrected positions of the MS-DArT Tags
input1  input1  bed file with all MS-DArT Tags
marks_with_msp_bigger_than_0.R -h | --help  show this message
" -> doc

# load the docopt library
require(docopt)
require(gdata)
require(dplyr)
set.seed(1407)
# retrieve the command-line arguments
opts <- docopt(doc)

# save.image("marks_with_msp_bigger_than_0.rda")
# load("marks_with_msp_bigger_than_0.rda")
dados <- read.table(opts$`<input1>`,
                    header = T,
                    sep = "\t",
                    check.names = F, row.names = 1)

# Filter all sites with counts larger than one in the samples coming from a methylation-insensitive library.

## first, select columns from the methylation insensitive library
samples_from_ms_insenstive <- strsplit(opts$mspI, split = ",")[[1]]
samples_from_ms_insenstive <- paste0("mapping/", samples_from_ms_insenstive, "_combined.bam")

dados <- dados[, colnames(dados) %in% samples_from_ms_insenstive ]
dados <- cbind(dados[,1], dados)

## Then, filter the columns
dados2 <-dados %>%
    filter_all(all_vars(.>0))

dados3 <- as.data.frame(row.names(dados2))

write.table(dados3,
            paste(opts$O1),
            row.names = F,
            quote = F,
            sep = "\t")
