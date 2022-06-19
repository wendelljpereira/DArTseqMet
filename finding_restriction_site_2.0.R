"Usage: finding_restriction_site_2.0.R (--out1 <O1>) <input1> <input2>
-h --help    show this help
--out1  name1   bed file with the position of all mspI restriction sites in the genome
<input1>  input1  fasta file with the genome
<input2>  input2  fasta file with the restriction site sequence of the enzymes
finding_restriction_site_2.0.R -h | --help  show this message
" -> doc

# load the docopt library
suppressMessages ( require("Biostrings") )
suppressMessages ( require("docopt") )

# retrieve the command-line arguments
opts <- docopt::docopt(doc)

save.image("test.rda")

genome <- Biostrings::readDNAStringSet(opts$`<input1>`)
restriction_site_list <- Biostrings::readDNAStringSet(opts$`<input2>`)

# Loads the functions 
find_restriction_sites <- function(i,
                                   restriction_site_list=restriction_site_list,
                                   target_sequence=target_sequence, # genome[i][[1]]
                                   target_sequence_name = target_sequence_name # names(genome[1])
) {
    
    enzime_name <- names(restriction_site_list)[i]
    restriction_site <- restriction_site_list[[i]]
    
    # search in plus strand
    plus_matches <- matchPattern(restriction_site, target_sequence)
    rep.int(enzime_name, length(plus_matches))
    
    data.frame(rep(target_sequence_name, length(plus_matches)),
               start(plus_matches) - 1,
               end(plus_matches),
               rep(enzime_name, length(plus_matches)),
               rep(0, length(plus_matches)),
               rep("+", length(plus_matches)))
    
}

# For each enzime, look for the restrictions sites in all chromosomes
all_stances_of_the_site <- data.frame()
for( RE in 1:length( names(restriction_site_list) ) ) {
    
    ## For each chromosome, find the instances where the restriction site is.
    for (chr in 1:length( names(genome) ) ) {
        #@@@     for ( chr in 1:2 ) {    
        RE_sites <- find_restriction_sites(i = RE,
                                           restriction_site_list=restriction_site_list,
                                           target_sequence = genome[chr][[1]], # genome[i][[1]]
                                           target_sequence_name = names(genome[chr]) # names(genome[1]) 
        )
        
        all_stances_of_the_site <- rbind(all_stances_of_the_site, RE_sites)
        
    }
    
}

sites <- all_stances_of_the_site
colnames(sites) <- c("V1", "V2", "V3", "V4", "V5", "V6")

# Writes the bed file containing the restriction sites.
sites_minus <- sites
sites_minus$V6 <- "-"

sites <- rbind(sites, sites_minus)

# Ordena o arquivo pela posicao
sites <- sites[with(sites, order(V1, V2)), ]

write.table(sites,
            paste(opts$O1),
            col.names = F,
            row.names = F,
            quote = F,
            sep = "\t")
