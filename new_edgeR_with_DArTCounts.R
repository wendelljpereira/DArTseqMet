require(ggplot2)
require(edgeR)
require(plyr)
require(gdata)
require(tidyr)
set.seed(1407)

# save.image("edgeR_with_DArTCounts.rda")
# load("edgeR_with_DArTCounts.rda")

fdr = as.numeric( snakemake@params[["fdr"]] )
log_fc = as.numeric( snakemake@params[["log_fold_change"]] )
prefix = snakemake@params[["prefix"]]
number_of_tec_rep <- as.numeric( snakemake@params[["number_of_tec_rep"]] )
min_msp = as.numeric(snakemake@params[["min_msp"]])

data <- read.table(as.character(snakemake@input[2]),
                   header = T,
                   sep = ",",
                   check.names = F)

marcas <- as.character(data[, 1]) # Selects the names of the sites
prefixo <- as.character(prefix)

## Define the names on each group
Grouping_var <- as.data.frame( strsplit(snakemake@params[["Grouping_var"]], split = ",")[[1]] )
Grouping_var <- tidyr::separate(Grouping_var, col = 1, into = c("sample", "group"), sep = "-" )

groups_names <- unique(Grouping_var$group)

# For each listed groups, find the DNA methylations
for( gn in groups_names ) {
    
    sample_names <- Grouping_var[Grouping_var$group == gn, ]
    sample_names <- paste0("mapping/", unique(sample_names$sample), "_combined.bam")
    
    # Selects the name of the sample and tissue that should be used
    clone <- data[, colnames(data) %in% sample_names ]
    clone <- as.data.frame(clone, row.names = marcas)
    
    intersect <- read.table(paste(snakemake@input[1]),
                            header = T,
                            quote = "\"",
                            sep = "\t",
                            check.names = F,
                            colClasses = "character")
    
    intersect <- as.character(intersect[, 1])
    
    clone <- clone[rownames(clone) %in% intersect, ]
    
    # Filters by the minimum counts to be considerated as a true site
    filtred_clone <- data.frame()
    
    for (i in 1:nrow(clone) ) {
        
        count_average <- clone[i, ]
        if ( (sum(count_average[1, ]) / length(count_average)) >= min_msp ) {
            
            filtred_clone <- rbind(filtred_clone, clone[i, ])
            
        }
    }
    
    msp_10_sites <- as.data.frame( rownames(filtred_clone) )
    clone <- filtred_clone
    
    samples_from_ms_insenstive <- strsplit(snakemake@params[["samples_from_ms_insenstive"]], split = ",")[[1]]
    
    groups <- data.frame( "samples" = colnames(clone) )
    groups <- groups %>%
        tidyr::separate(col = samples, into = c("a", "sample"), sep = "\\/") %>%
        tidyr::separate(col = sample, into = c("sample", "b"), sep = "\\_") %>%
        dplyr::select(sample)
    
    groups$group <- ifelse(groups$sample %in% samples_from_ms_insenstive, "ms", "hp" )
    
    # Loads the processed data to edgeR object
    data_edgeR <- DGEList(counts = clone, group = groups$group)
    
    # Determines the dispersion and save the values in a new file
    data_edgeR <- estimateDisp(data_edgeR)
    
    # Determines the significantly different counts between MspI and HpaII
    tags <- exactTest(data_edgeR, pair = c("hp", "ms"))
    
    tTags <- topTags(tags, n = NULL)
    results <- as.data.frame(tTags)
    
    results_sig <- subset(results, results$FDR <= fdr)
    results_sig <- subset(results_sig, results_sig$logFC >= log_fc)
    
    if ( nrow( results_sig ) == 0 ) {
        
        print("There is no differentialy expressed sites for this sample.")
        
    } else {
        
        marcks_sig <- rownames(results_sig)
        Data_DE <- as.data.frame(marcks_sig)
        
        colnames(Data_DE) <- gn
        
        if ( file.exists( paste( prefixo, "DE_marks.txt", sep = "_" ) ) == "FALSE" ) {
            
            write.table(Data_DE,
                        file = paste(prefixo, "DE_marks.txt", sep = "_"),
                        sep = "\t",
                        quote = F,
                        row.names = F,
                        col.names = T)
            
        } else if ( file.exists( paste( prefixo, "DE_marks.txt", sep = "_" ) ) == "TRUE" ) {
            
            Data_DE_g <- read.table(paste(prefixo, "DE_marks.txt", sep = "_"),
                                    header = T,
                                    sep = "\t")
            
            Data_DE_g <- cbindX(Data_DE_g, Data_DE)
            
            write.table(Data_DE_g,
                        file = paste(prefixo, "DE_marks.txt", sep = "_"),
                        sep = "\t",
                        quote = F,
                        row.names = F,
                        col.names = T)
            
        }
    }
}
