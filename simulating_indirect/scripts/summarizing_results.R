#log <- file(snakemake@log[[1]], open="wt")
#sink(log)
#sink(log, type="message")

##code to determine the number of TP, FP, TN, FN

#From the paper: egene is a TP if simulate eQTL is among the significant eQTLs
#Note that the FN includes eGenes with no significant eSNPs (type 1) 
#and eGenes with significant eSNPs, but for which the correct eSNP was not 
#significant (type 2).
#Each simulated gene not assigned an eQTL effect was considered a 
#true negative (TN) if it had no significant eSNPs and a 
#false positive (FP) if it had a significant hit.

#####Libraries#####
library(optparse)
library(qvalue)
library(tidyverse)
library(data.table)

#####Variable loading#####
option_list <- list(
  make_option("--key", type="character"),
  make_option("--results", type="character"),
  make_option("--effect_size", type="double"),
  make_option("--samp_n", type="integer"),
  make_option("--cell_count", type="integer"),
  make_option("--rep", type="integer"),
  make_option("--rand", type="integer")
)

opt <- parse_args(OptionParser(option_list=option_list))

outputfile <- paste0("results_celltype_A/", as.character(opt$effect_size), "_sample",
  as.character(opt$samp_n), "_cell", as.character(opt$cell_count), "_rep",
  as.character(opt$rep), "_rand", as.character(opt$rand), "/qtl_stats.txt")

#first going to bring the egenes into positive and negative
key_name <- opt$key
results_name <- opt$results

print("start")
#going to go through the key

#first, going to add a column with the adjusted p value to the qtl result table
#filter out all rows that don't have pvalue below the threshold
#get unique of the genes to see all the genes with simulated 
key <- readRDS(key_name)
key <- key$key
key$geneID <- rownames(key)
key <- key[
  !is.na(key$eSNP.ID) &
    key$eQTL.group %in% c("global", "Group1"),   # global eQTLs + Group1-specific ones
  c("geneID", "eSNP.ID", "eQTL.EffectSize")
]

results <- data.table::fread(results_name)
results$qval <- qvalue(results$p_value)$qvalue
sig_results <- results[results$qval <= 0.05, ]
sig_egenes <- unique(sig_results$feature_id)

getter <- function(gene_name, esnp) {
  if (!is.na(esnp)) {
    if (!(any(results$feature_id == gene_name & 
              results$snp_id == esnp, na.rm = TRUE))) {
      return("Not tested")
    }
  }
  
  if (gene_name %in% sig_egenes) { #seeing if the gene is signficant
    #so subsetting all the snps in the egene of choice
    sig_snps <- sig_results[sig_results$feature_id == gene_name]
    if (!is.na(esnp) &&  esnp %in% sig_snps$snp_id) { #seeing that there is the right eQTL in it
      return("TP")
    }  
    else {return("FP")}
  }
  else { 
    if (is.na(esnp)){ #i.e. model didn't say gene is significant and that is right 
      return("TN")
    }
    else {return("FN")}
  }
}

stats_results <- mapply(getter, rownames(key), key$eSNP.ID)


#getting beta correlation values
key$gene_name <- rownames(key)
key_sub <- key[!(is.na(key$'eSNP.ID')), c("gene_name", "eSNP.ID", "eQTL.EffectSize")]
result_sub <- results[,c("feature_id", "snp_id", "beta")]

merged_df <- merge(key_sub, 
      result_sub, 
      by.x = c("eSNP.ID", "gene_name"), 
      by.y = c("snp_id", "feature_id"))

beta_corr <- cor(abs(merged_df$eQTL.EffectSize), 
                 abs(merged_df$beta))

#output files 
output_table <- data.frame("effect_size" = opt$effect_size,
           "sample_size" = opt$samp_n,
           "cell_count" = opt$cell_count,
           "rep" = opt$rep,
           "FN" = sum(stats_results == "FN"),
           "FP" = sum(stats_results == "FP"),
           "TN" = sum(stats_results == "TN"),
           "TP" = sum(stats_results == "TP"),
           "beta_corr" = beta_corr)

write.table(output_table, file = outputfile, sep = ",", 
            row.names = FALSE, col.names = TRUE,
            quote = FALSE)

